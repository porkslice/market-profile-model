// z_Bayes_Bridge.cpp  — Sierra ↔ Python bridge + baseline profile visual
#include "sierrachart.h"
#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cstdio>

SCDLLName("Bayesian Profile – Bridge")

// --------- small helpers ---------
static std::string JoinPath(const SCString& a, const char* b){
  SCString out = a; if(out.GetLength() && out[out.GetLength()-1] != '\\') out += "\\";
  out += b; return std::string(out.GetChars());
}
static bool AtomicWrite(const std::string& path, const std::string& content){
  std::string tmp = path + ".tmp";
  FILE* f = nullptr;
#ifdef _MSC_VER
  fopen_s(&f, tmp.c_str(), "wb");
#else
  f = fopen(tmp.c_str(), "wb");
#endif
  if(!f) return false;
  fwrite(content.data(), 1, content.size(), f);
  fclose(f);
  // Replace target
  remove(path.c_str());
  return rename(tmp.c_str(), path.c_str()) == 0;
}
static SCDateTime NowDT(const SCStudyInterfaceRef& sc){ return sc.BaseDateTimeIn[sc.ArraySize-1]; }

// --------- Study ---------
SCSFExport scsf_BayesProfileBridge(SCStudyInterfaceRef sc)
{
  // Inputs
  SCInputRef In_GridTicks       = sc.Input[0];
  SCInputRef In_SigmaBins       = sc.Input[1];
  SCInputRef In_SigmaMin        = sc.Input[2];
  SCInputRef In_SigmaMax        = sc.Input[3];
  SCInputRef In_VAPct           = sc.Input[4];
  SCInputRef In_OR_Min          = sc.Input[5];
  SCInputRef In_RefreshSec      = sc.Input[6];
  SCInputRef In_DrawPriorRefs   = sc.Input[7];
  SCInputRef In_OutDir          = sc.Input[8];

  // Subgraphs (visuals)
  SCSubgraphRef SG_POC          = sc.Subgraph[0];
  SCSubgraphRef SG_VA_Up        = sc.Subgraph[1];
  SCSubgraphRef SG_VA_Dn        = sc.Subgraph[2];
  SCSubgraphRef SG_PostPOC      = sc.Subgraph[3]; // posterior POC (Python)
  SCSubgraphRef SG_PostVA_Up    = sc.Subgraph[4];
  SCSubgraphRef SG_PostVA_Dn    = sc.Subgraph[5];
  SCSubgraphRef SG_LamMax       = sc.Subgraph[6]; // λmax sparkline (posterior context)

  // Persist
  int& p_SessionStartIdx   = sc.GetPersistentInt(1);
  int& p_OR_EndIdx         = sc.GetPersistentInt(2);
  int& p_LastRefreshEpoch  = sc.GetPersistentInt(3);
  int& p_DayKey            = sc.GetPersistentInt(4);
  float& p_SigmaActive     = sc.GetPersistentFloat(1);

  // Prior references (persist)
  double& Prior_VAH = sc.GetPersistentDouble(1);
  double& Prior_VAL = sc.GetPersistentDouble(2);
  double& Prior_POC = sc.GetPersistentDouble(3);

  if (sc.SetDefaults)
  {
    sc.GraphName = "Bayesian Profile – Bridge";
    sc.AutoLoop = 0;
    sc.GraphRegion = 0;
    sc.MaintainVolumeAtPriceData = 1;
    sc.UpdateAlways = 1;

    In_GridTicks.Name = "Profile Grid Size (ticks)";  In_GridTicks.SetInt(2); In_GridTicks.SetIntLimits(1,5);
    In_SigmaBins.Name = "Smoothing Sigma (bins)";     In_SigmaBins.SetFloat(2.0f);
    In_SigmaMin.Name  = "Sigma Min (bins)";           In_SigmaMin.SetFloat(1.0f);
    In_SigmaMax.Name  = "Sigma Max (bins)";           In_SigmaMax.SetFloat(4.0f);
    In_VAPct.Name     = "Value Area Percent (0..1)";  In_VAPct.SetFloat(0.70f);
    In_OR_Min.Name    = "Opening Range Minutes";      In_OR_Min.SetInt(5); In_OR_Min.SetIntLimits(1,30);
    In_RefreshSec.Name= "Refresh Seconds";            In_RefreshSec.SetInt(10); In_RefreshSec.SetIntLimits(2,60);
    In_DrawPriorRefs.Name = "Draw Prior VAH/VAL/POC"; In_DrawPriorRefs.SetYesNo(true);
    In_OutDir.Name    = "Output Subfolder (under Data Files Folder)";
    In_OutDir.SetString("bayes");

    SG_POC.Name="POC (fallback)"; SG_POC.DrawStyle = DRAWSTYLE_LINE; SG_POC.PrimaryColor = RGB(220,220,220); SG_POC.LineWidth=2;
    SG_VA_Up.Name="VA Upper (fallback)"; SG_VA_Up.DrawStyle=DRAWSTYLE_LINE; SG_VA_Up.PrimaryColor=RGB(0,200,200); SG_VA_Up.LineStyle=LINESTYLE_DASH;
    SG_VA_Dn.Name="VA Lower (fallback)"; SG_VA_Dn.DrawStyle=DRAWSTYLE_LINE; SG_VA_Dn.PrimaryColor=RGB(0,200,200); SG_VA_Dn.LineStyle=LINESTYLE_DASH;

    SG_PostPOC.Name="POC (posterior)"; SG_PostPOC.DrawStyle=DRAWSTYLE_LINE; SG_PostPOC.PrimaryColor=RGB(255,255,255); SG_PostPOC.LineWidth=2;
    SG_PostVA_Up.Name="VA Upper (posterior)"; SG_PostVA_Up.DrawStyle=DRAWSTYLE_LINE; SG_PostVA_Up.PrimaryColor=RGB(255,255,255); SG_PostVA_Up.LineStyle=LINESTYLE_DOT;
    SG_PostVA_Dn.Name="VA Lower (posterior)"; SG_PostVA_Dn.DrawStyle=DRAWSTYLE_LINE; SG_PostVA_Dn.PrimaryColor=RGB(255,255,255); SG_PostVA_Dn.LineStyle=LINESTYLE_DOT;

    SG_LamMax.Name="λmax"; SG_LamMax.DrawStyle=DRAWSTYLE_LINE; SG_LamMax.PrimaryColor=RGB(255,180,0);

    return;
  }

  const int last = sc.ArraySize - 1;
  if (last < 10) return;

  const int todayKey = sc.GetTradingDayDate(sc.BaseDateTimeIn[last]);
  if (p_DayKey != todayKey || sc.IsNewTradingDay(sc.BaseDateTimeIn[last]))
  {
    // New RTH session
    p_DayKey = todayKey;
    p_SessionStartIdx = last;
    for (int i = last; i >= 0; --i){
      if (sc.GetTradingDayDate(sc.BaseDateTimeIn[i]) != todayKey) break;
      p_SessionStartIdx = i;
    }
    SCDateTime orEnd = sc.BaseDateTimeIn[p_SessionStartIdx] + SCDateTime::MINUTES(In_OR_Min.GetInt());
    p_OR_EndIdx = sc.GetContainingIndexForSCDateTime(sc.ChartNumber, orEnd);
    p_SigmaActive = (float)In_SigmaBins.GetFloat();

    // Build prior-day references (simple VAP)
    Prior_VAH = Prior_VAL = Prior_POC = 0.0;
    if (In_DrawPriorRefs.GetYesNo())
    {
      int prevEnd=-1, prevStart=-1;
      for(int i=p_SessionStartIdx-1; i>=0; --i){ if(prevEnd==-1) prevEnd=i; if(sc.GetTradingDayDate(sc.BaseDateTimeIn[i]) != todayKey){ prevStart=i+1; break; } }
      if(prevStart>=0 && prevEnd>prevStart){
        std::map<int64_t,double> vap;
        for(int bi=prevStart; bi<=prevEnd; ++bi){
          const int cnt = sc.VolumeAtPriceForBars->GetSizeAtBarIndex(bi);
          s_VolumeAtPriceV2* e=nullptr;
          for(int vi=0; vi<cnt; ++vi){ sc.VolumeAtPriceForBars->GetVAPElementAtIndex(bi,vi,&e); if(!e) continue; vap[e->PriceInTicks] += (double)e->Volume; }
        }
        if(!vap.empty()){
          std::vector<int64_t> px; std::vector<double> vv; px.reserve(vap.size()); vv.reserve(vap.size());
          double total=0; for(auto& kv:vap){ px.push_back(kv.first); vv.push_back(kv.second); total+=kv.second; }
          size_t pocI = (size_t)(std::max_element(vv.begin(),vv.end())-vv.begin());
          Prior_POC = sc.TickSize * (double)px[pocI];
          // VA by POC expansion
          double target = total * In_VAPct.GetFloat(); size_t L=pocI, R=pocI; double acc=vv[pocI];
          while(acc<target && (L>0 || R+1<vv.size())){
            double left=(L>0)?vv[L-1]:-1.0, right=(R+1<vv.size())?vv[R+1]:-1.0;
            if(right>=left && R+1<vv.size()){ R++; acc+=right; }
            else if(L>0){ L--; acc+=left; }
            else break;
          }
          Prior_VAL = sc.TickSize * (double)px[L];
          Prior_VAH = sc.TickSize * (double)px[R];

          // draw lines
          s_UseTool t; t.Clear(); t.ChartNumber=sc.ChartNumber; t.AddMethod=UTAM_ADD_OR_ADJUST;
          t.Region=sc.GraphRegion; t.LineStyle=LINESTYLE_DOT; t.LineWidth=1;
          t.Color=RGB(255,0,200); t.DrawingType=DRAWING_HORIZONTALLINE; t.BeginValue=Prior_VAH; sc.UseTool(t);
          t.BeginValue=Prior_VAL; sc.UseTool(t);
          t.Color=RGB(200,200,200); t.LineStyle=LINESTYLE_DASH; t.BeginValue=Prior_POC; sc.UseTool(t);
        }
      }
    }
  }

  // Throttle
  const int refresh = In_RefreshSec.GetInt();
  const int nowSec = (int)sc.BaseDateTimeIn[last].GetTimeInSeconds();
  if (nowSec - p_LastRefreshEpoch < refresh && sc.UpdateStartIndex!=0) return;
  p_LastRefreshEpoch = nowSec;

  // ----- Build session VAP (grid) -----
  const int grid = In_GridTicks.GetInt();
  std::map<int64_t,double> vap; vap.clear();
  for(int bi=p_SessionStartIdx; bi<=last; ++bi){
    const int cnt = sc.VolumeAtPriceForBars->GetSizeAtBarIndex(bi);
    s_VolumeAtPriceV2* e=nullptr;
    for(int vi=0; vi<cnt; ++vi){ sc.VolumeAtPriceForBars->GetVAPElementAtIndex(bi,vi,&e); if(!e) continue;
      int64_t g = (e->PriceInTicks / grid) * (int64_t)grid;
      vap[g] += (double)e->Volume;
    }
  }
  if(vap.size()<5) return;
  std::vector<int64_t> px; std::vector<double> vv; px.reserve(vap.size()); vv.reserve(vap.size());
  for(auto& kv:vap){ px.push_back(kv.first); vv.push_back(kv.second); }

  // ----- Smooth (fallback) -----
  auto smooth = [](const std::vector<double>& x, double s){
    if(x.empty()) return x;
    int n=(int)x.size(); std::vector<double> y(n,0.0);
    int R=(int)std::ceil(3*s); double denom=2*s*s;
    for(int i=0;i<n;++i){ double acc=0, wsum=0;
      int L=std::max(0,i-R), U=std::min(n-1,i+R);
      for(int j=L;j<=U;++j){ double d=j-i; double w=exp(-(d*d)/denom); acc+=w*x[j]; wsum+=w; }
      y[i]=(wsum>0?acc/wsum:x[i]);
    }
    return y;
  };
  float sMin=(float)In_SigmaMin.GetFloat(), sMax=(float)In_SigmaMax.GetFloat();
  float sNew=(float)In_SigmaBins.GetFloat(); sNew=std::min(std::max(sNew,sMin),sMax);
  if(p_SigmaActive<=0) p_SigmaActive=sNew; else{
    float rel=fabsf(sNew-p_SigmaActive)/std::max(1e-6f,p_SigmaActive);
    if(rel>=0.10f) p_SigmaActive=sNew;
  }
  std::vector<double> sm = smooth(vv, p_SigmaActive);
  double total=0; for(double v:sm) total+=std::max(0.0,v);
  if(total<=0) return;
  for(double& v:sm) v = std::max(1e-20, v/total);

  // POC & VA (fallback)
  size_t pocI = (size_t)(std::max_element(sm.begin(),sm.end())-sm.begin());
  double POC_fallback = sc.TickSize*(double)px[pocI];
  double target = std::min(0.99,std::max(0.50,(double)In_VAPct.GetFloat()));
  size_t L=pocI,R=pocI; double acc=sm[pocI];
  while(acc<target && (L>0 || R+1<sm.size())){
    double left=(L>0)?sm[L-1]:-1.0, right=(R+1<sm.size())?sm[R+1]:-1.0;
    if(right>=left && R+1<sm.size()){R++; acc+=std::max(0.0,right);}
    else if(L>0){L--; acc+=std::max(0.0,left);}
    else break;
  }
  double VAlo_fallback = sc.TickSize*(double)px[L];
  double VAhi_fallback = sc.TickSize*(double)px[R];

  for(int i=sc.UpdateStartIndex;i<=last;++i){ SG_POC[i]=POC_fallback; SG_VA_Dn[i]=VAlo_fallback; SG_VA_Up[i]=VAhi_fallback; }

  // ----- Draw OR box -----
  int orStart=p_SessionStartIdx, orEnd=std::min(p_OR_EndIdx,last);
  if(orEnd>orStart){
    double orHi=sc.High[orStart], orLo=sc.Low[orStart];
    for(int i=orStart;i<=orEnd;++i){ orHi=std::max(orHi,sc.High[i]); orLo=std::min(orLo,sc.Low[i]); }
    s_UseTool t; t.Clear(); t.ChartNumber=sc.ChartNumber; t.AddMethod=UTAM_ADD_OR_ADJUST;
    t.DrawingType=DRAWING_RECTANGLEHIGHLIGHT; t.Color=RGB(70,130,180); t.TransparencyLevel=85;
    t.BeginIndex=orStart; t.EndIndex=orEnd; t.BeginValue=orLo; t.EndValue=orHi; t.Region=sc.GraphRegion; sc.UseTool(t);
  }

  // ----- WRITE FEEDS -----
  const std::string base = JoinPath(sc.DataFilesFolder(), In_OutDir.GetString());
  // Ensure folder exists (silent ok)
  CreateDirectoryA(base.c_str(), NULL);

  // feed_bins.csv: timestamp,grid_ticks,tick_size,price_ticks,volume
  {
    std::ostringstream oss;
    SCDateTime now = NowDT(sc);
    int date = now.GetDate(); int time = now.GetTimeInSeconds();
    oss << "timestamp,date,time,grid_ticks,tick_size,price_ticks,volume\n";
    for(size_t i=0;i<px.size();++i){
      oss << (double)now << "," << date << "," << time << "," << grid << "," << sc.TickSize
          << "," << (long long)px[i] << "," << (long long)vv[i] << "\n";
    }
    AtomicWrite(base + "\\feed_bins.csv", oss.str());
  }

  // feed_features.csv: one row snapshot with session refs + OR + prior refs + sigma
  {
    SCDateTime now = NowDT(sc);
    int date = now.GetDate(); int time = now.GetTimeInSeconds();
    double openPx = sc.Open[p_SessionStartIdx];
    // Simple ADR20 proxy: average true range of last 20 RTH days (fallback: use prior range only if needed)
    double ADR20 = 0.0; int days=0, idx = p_SessionStartIdx-1;
    while(idx>0 && days<20){
      int d = sc.GetTradingDayDate(sc.BaseDateTimeIn[idx]);
      int j=idx; double hi=sc.High[idx], lo=sc.Low[idx];
      while(j>0 && sc.GetTradingDayDate(sc.BaseDateTimeIn[j])==d){ hi=std::max(hi,sc.High[j]); lo=std::min(lo,sc.Low[j]); --j; }
      ADR20 += (hi-lo); days++; idx=j;
    }
    if(days>0) ADR20/=days;

    double prevHi = (p_SessionStartIdx>0)? sc.High[p_SessionStartIdx-1]:0.0;
    double prevLo = (p_SessionStartIdx>0)? sc.Low[p_SessionStartIdx-1]:0.0;
    double gap = (ADR20>0? (openPx - sc.Close[p_SessionStartIdx-1]) / ADR20 : 0.0);

    // OR bounds
    double orHi=sc.High[orStart], orLo=sc.Low[orStart];
    for(int i=orStart;i<=std::min(orEnd,last);++i){ orHi=std::max(orHi,sc.High[i]); orLo=std::min(orLo,sc.Low[i]); }

    std::ostringstream oss;
    oss << "timestamp,date,time,grid_ticks,tick_size,session_open,prior_poc,prior_val,prior_vah,prev_high,prev_low,adr20,"
           "or_low,or_high,sigma_bins,profile_poc_fallback,va_low_fallback,va_high_fallback\n";
    oss << (double)now << "," << date << "," << time << "," << grid << "," << sc.TickSize << ","
        << openPx << "," << Prior_POC << "," << Prior_VAL << "," << Prior_VAH << ","
        << prevHi << "," << prevLo << "," << ADR20 << ","
        << orLo << "," << orHi << "," << p_SigmaActive << ","
        << POC_fallback << "," << VAlo_fallback << "," << VAhi_fallback << "\n";
    AtomicWrite(base + "\\feed_features.csv", oss.str());
  }

  // ----- READ RESULTS -----
  {
    std::string path = base + "\\results.csv";
    std::ifstream fin(path.c_str());
    if(fin.good()){
      std::string line, lastline;
      std::getline(fin,line); // header
      while(std::getline(fin,line)) if(!line.empty()) lastline = line;
      fin.close();
      if(!lastline.empty()){
        // Parse minimal fields
        // header must be: timestamp,post_poc,post_va_low,post_va_high,post_poc_lo50,post_poc_hi50,post_va_low50,post_va_high50,lam_max,w1,w2,w3,w4
        std::stringstream ss(lastline);
        std::string tok; std::vector<std::string> cols;
        while(std::getline(ss,tok,',')) cols.push_back(tok);
        if(cols.size()>=10){
          double postPOC = atof(cols[1].c_str());
          double postVAL = atof(cols[2].c_str());
          double postVAH = atof(cols[3].c_str());
          double lammax  = atof(cols[8].c_str());
          for(int i=sc.UpdateStartIndex;i<=last;++i){
            SG_PostPOC[i]=postPOC; SG_PostVA_Dn[i]=postVAL; SG_PostVA_Up[i]=postVAH;
            SG_LamMax[i]=lammax;
          }
          // Lightweight HUD
          s_UseTool text; text.Clear(); text.ChartNumber=sc.ChartNumber; text.AddMethod=UTAM_ADD_OR_ADJUST;
          text.DrawingType=DRAWING_TEXT; text.FontSize=10; text.Color=RGB(255,255,255); text.Region=sc.GraphRegion;
          SCString hud; hud.Format("Posterior: POC %.2f | VA [%.2f, %.2f] | λmax %.2f", postPOC, postVAL, postVAH, lammax);
          text.Text = hud; text.BeginIndex=last; text.BeginValue=sc.LastTradePrice; sc.UseTool(text);
        }
      }
    }
  }
}
