// sierra_chart_marketprofilemodel.cpp — Bridge + priors + bands + flow markers + weights HUD
#include "sierrachart.h"
#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cstdio>
#include <cmath>

SCDLLName("Bayesian Market-Profile – Bridge+Viz")

// ---------- helpers ----------
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
  remove(path.c_str());
  return rename(tmp.c_str(), path.c_str()) == 0;
}
static SCDateTime NowDT(const SCStudyInterfaceRef& sc){ return sc.BaseDateTimeIn[sc.ArraySize-1]; }

// ---------- study ----------
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
  SCInputRef In_FlowLookback    = sc.Input[9];
  SCInputRef In_ImbalanceThresh = sc.Input[10];
  SCInputRef In_DrawBands       = sc.Input[11];

  // Subgraphs
  SCSubgraphRef SG_POC          = sc.Subgraph[0];
  SCSubgraphRef SG_VA_Up        = sc.Subgraph[1];
  SCSubgraphRef SG_VA_Dn        = sc.Subgraph[2];
  SCSubgraphRef SG_PostPOC      = sc.Subgraph[3];
  SCSubgraphRef SG_PostVA_Up    = sc.Subgraph[4];
  SCSubgraphRef SG_PostVA_Dn    = sc.Subgraph[5];
  SCSubgraphRef SG_LamMax       = sc.Subgraph[6];

  // Persist
  int&   p_SessionStartIdx  = sc.GetPersistentInt(1);
  int&   p_OR_EndIdx        = sc.GetPersistentInt(2);
  int&   p_LastRefreshEpoch = sc.GetPersistentInt(3);
  int&   p_DayKey           = sc.GetPersistentInt(4);
  float& p_SigmaActive      = sc.GetPersistentFloat(1);

  // Prior refs
  double& Prior_VAH = sc.GetPersistentDouble(1);
  double& Prior_VAL = sc.GetPersistentDouble(2);
  double& Prior_POC = sc.GetPersistentDouble(3);

  // Opening type sticky
  int&   p_OpenTypeSticky   = sc.GetPersistentInt(5);   // 1=OD,2=OTD,3=ORR,4=OAIR,5=OAOR
  float& p_OpenTypeProb     = sc.GetPersistentFloat(2);
  float& p_LastOpenProb     = sc.GetPersistentFloat(3);
  int&   p_HadORBreak       = sc.GetPersistentInt(6);
  int&   p_ORFirstBreakIdx  = sc.GetPersistentInt(7);

  if (sc.SetDefaults)
  {
    sc.GraphName = "Bayesian Market-Profile – Bridge+Viz";
    sc.AutoLoop = 0; sc.GraphRegion = 0;
    sc.MaintainVolumeAtPriceData = 1; sc.UpdateAlways = 1;

    In_GridTicks.Name = "Profile Grid Size (ticks)";  In_GridTicks.SetInt(2); In_GridTicks.SetIntLimits(1,5);
    In_SigmaBins.Name = "Smoothing Sigma (bins)";     In_SigmaBins.SetFloat(2.0f);
    In_SigmaMin.Name  = "Sigma Min (bins)";           In_SigmaMin.SetFloat(1.0f);
    In_SigmaMax.Name  = "Sigma Max (bins)";           In_SigmaMax.SetFloat(4.0f);
    In_VAPct.Name     = "Value Area Percent (0..1)";  In_VAPct.SetFloat(0.70f);
    In_OR_Min.Name    = "Opening Range Minutes";      In_OR_Min.SetInt(5); In_OR_Min.SetIntLimits(1,30);
    In_RefreshSec.Name= "Refresh Seconds";            In_RefreshSec.SetInt(10); In_RefreshSec.SetIntLimits(2,60);
    In_DrawPriorRefs.Name = "Draw Prior VAH/VAL/POC"; In_DrawPriorRefs.SetYesNo(true);
    In_OutDir.Name    = "Output Subfolder (under Data Files Folder)"; In_OutDir.SetString("bayes");
    In_FlowLookback.Name = "Flow Lookback Bars at Levels"; In_FlowLookback.SetInt(10); In_FlowLookback.SetIntLimits(3,50);
    In_ImbalanceThresh.Name = "Flow Imbalance Threshold (ratio)"; In_ImbalanceThresh.SetFloat(1.8f);
    In_DrawBands.Name = "Draw VA 50% Credible Band"; In_DrawBands.SetYesNo(true);

    SG_POC.Name="POC (fallback)"; SG_POC.DrawStyle=DRAWSTYLE_LINE; SG_POC.PrimaryColor=RGB(200,200,200); SG_POC.LineWidth=2;
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

  // New session?
  const int todayKey = sc.GetTradingDayDate(sc.BaseDateTimeIn[last]);
  bool is_new_td = false;
  if (p_DayKey != todayKey){
    is_new_td = true;
  } else if (last > 0) {
    int prevKey = sc.GetTradingDayDate(sc.BaseDateTimeIn[last-1]);
    is_new_td = (prevKey != todayKey);
  }
  if (is_new_td)
  {
    p_DayKey = todayKey;
    p_SessionStartIdx = last;
    for (int i = last; i >= 0; --i){
      if (sc.GetTradingDayDate(sc.BaseDateTimeIn[i]) != todayKey) break;
      p_SessionStartIdx = i;
    }
    SCDateTime orEnd = sc.BaseDateTimeIn[p_SessionStartIdx] + SCDateTime::MINUTES(In_OR_Min.GetInt());
    p_OR_EndIdx = sc.GetContainingIndexForSCDateTime(sc.ChartNumber, orEnd);
    p_SigmaActive = (float)In_SigmaBins.GetFloat();

    // reset opening type
    p_OpenTypeSticky = 0; p_OpenTypeProb = 0.0f; p_LastOpenProb = 0.0f;
    p_HadORBreak = 0; p_ORFirstBreakIdx = -1;

    // Prior-day refs
    Prior_VAH = Prior_VAL = Prior_POC = 0.0;
    if (In_DrawPriorRefs.GetYesNo())
    {
      int prevEnd=-1, prevStart=-1;
      for(int i=p_SessionStartIdx-1; i>=0; --i){ if(prevEnd==-1) prevEnd=i; if(sc.GetTradingDayDate(sc.BaseDateTimeIn[i]) != todayKey){ prevStart=i+1; break; } }
      if(prevStart>=0 && prevEnd>prevStart){
        std::map<int64_t,double> vap;
        for(int bi=prevStart; bi<=prevEnd; ++bi){
          int cnt = sc.VolumeAtPriceForBars->GetSizeAtBarIndex(bi);
          s_VolumeAtPriceV2* e=nullptr;
          for(int vi=0; vi<cnt; ++vi){ sc.VolumeAtPriceForBars->GetVAPElementAtIndex(bi,vi,&e); if(!e) continue; vap[e->PriceInTicks] += (double)e->Volume; }
        }
        if(!vap.empty()){
          std::vector<int64_t> px; std::vector<double> vv; px.reserve(vap.size()); vv.reserve(vap.size());
          double total=0; for(auto& kv:vap){ px.push_back(kv.first); vv.push_back(kv.second); total+=kv.second; }
          size_t pocI = (size_t)(std::max_element(vv.begin(),vv.end())-vv.begin());
          Prior_POC = sc.TickSize * (double)px[pocI];
          double target = total * In_VAPct.GetFloat(); size_t L=pocI, R=pocI; double acc=vv[pocI];
          while(acc<target && (L>0 || R+1<vv.size())){
            double left=(L>0)?vv[L-1]:-1.0, right=(R+1<vv.size())?vv[R+1]:-1.0;
            if(right>=left && R+1<vv.size()){ R++; acc+=right; }
            else if(L>0){ L--; acc+=left; }
            else break;
          }
          Prior_VAL = sc.TickSize*(double)px[L];
          Prior_VAH = sc.TickSize*(double)px[R];

          s_UseTool t; t.Clear(); t.ChartNumber=sc.ChartNumber; t.AddMethod=UTAM_ADD_OR_ADJUST; t.Region=sc.GraphRegion;
          t.DrawingType=DRAWING_HORIZONTALLINE; t.LineStyle=LINESTYLE_DOT; t.LineWidth=1; t.Color=RGB(255,0,200); t.BeginValue=Prior_VAH; sc.UseTool(t);
          t.BeginValue=Prior_VAL; sc.UseTool(t);
          t.Color=RGB(200,200,200); t.LineStyle=LINESTYLE_DASH; t.BeginValue=Prior_POC; sc.UseTool(t);
        }
      }
    }
  }

  // Throttle
  const int refresh = In_RefreshSec.GetInt();
  const int nowSec = sc.SecondsSinceStartOfDay(sc.BaseDateTimeIn[last]);
  if (nowSec - p_LastRefreshEpoch < refresh && sc.UpdateStartIndex!=0) return;
  p_LastRefreshEpoch = nowSec;

  // Session VAP (grid)
  const int grid = In_GridTicks.GetInt();
  std::map<int64_t,double> vap; vap.clear();
  for(int bi=p_SessionStartIdx; bi<=last; ++bi){
    int cnt = sc.VolumeAtPriceForBars->GetSizeAtBarIndex(bi);
    s_VolumeAtPriceV2* e=nullptr;
    for(int vi=0; vi<cnt; ++vi){
      sc.VolumeAtPriceForBars->GetVAPElementAtIndex(bi,vi,&e);
      if(!e) continue;
      int64_t g = (e->PriceInTicks / grid) * (int64_t)grid;
      vap[g] += (double)e->Volume;
    }
  }
  if(vap.size()<5) return;
  std::vector<int64_t> px; std::vector<double> vv; px.reserve(vap.size()); vv.reserve(vap.size());
  for(auto& kv:vap){ px.push_back(kv.first); vv.push_back(kv.second); }

  // Gaussian smoothing (fallback)
  auto smooth = [](const std::vector<double>& x, double s){
    if(x.empty()) return x; int n=(int)x.size(); std::vector<double> y(n,0.0);
    int R=(int)std::ceil(3*s); double denom=2*s*s;
    for(int i=0;i<n;++i){ double acc=0, wsum=0; int L=std::max(0,i-R), U=std::min(n-1,i+R);
      for(int j=L;j<=U;++j){ double d=j-i; double w=exp(-(d*d)/denom); acc+=w*x[j]; wsum+=w; }
      y[i]=(wsum>0?acc/wsum:x[i]);
    } return y;
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

  // Opening Range box
  int orStart=p_SessionStartIdx, orEnd=std::min(p_OR_EndIdx,last);
  double orHi=0, orLo=0;
  if(orEnd>orStart){
    orHi=sc.High[orStart]; orLo=sc.Low[orStart];
    for(int i=orStart;i<=orEnd;++i){ orHi=std::max(orHi,sc.High[i]); orLo=std::min(orLo,sc.Low[i]); }
    s_UseTool t; t.Clear(); t.ChartNumber=sc.ChartNumber; t.AddMethod=UTAM_ADD_OR_ADJUST;
    t.DrawingType=DRAWING_RECTANGLEHIGHLIGHT; t.Color=RGB(70,130,180); t.TransparencyLevel=85;
    t.BeginIndex=orStart; t.EndIndex=orEnd; t.BeginValue=orLo; t.EndValue=orHi; t.Region=sc.GraphRegion; sc.UseTool(t);
  }

  // Opening type heuristics (sticky + probability)
  // detect first OR break
  if (!p_HadORBreak && last > orEnd){
    for(int i=orEnd+1;i<=last;++i){
      if(sc.High[i] > orHi || sc.Low[i] < orLo){ p_HadORBreak=1; p_ORFirstBreakIdx=i; break; }
    }
  }
  // break-and-hold mins
  double heldMins=0.0;
  if (p_HadORBreak && p_ORFirstBreakIdx>0){
    bool outsideNow = (sc.LastTradePrice > orHi || sc.LastTradePrice < orLo);
    if(outsideNow){
      int t1 = sc.SecondsSinceStartOfDay(sc.BaseDateTimeIn[p_ORFirstBreakIdx]);
      int t2 = sc.SecondsSinceStartOfDay(sc.BaseDateTimeIn[last]);
      heldMins = (t2 - t1) / 60.0;
    }
  }
  // vwap & rotations (first 30m)
  double vwap=0, pv=0, vol=0; for(int i=p_SessionStartIdx;i<=last;++i){ pv += sc.Volume[i]*sc.HLCAvg[i]; vol += sc.Volume[i]; }
  if(vol>0) vwap=pv/vol;
  SCDateTime halfHour = sc.BaseDateTimeIn[orStart] + SCDateTime::MINUTES(30);
  int endRotIdx = sc.GetContainingIndexForSCDateTime(sc.ChartNumber, halfHour); endRotIdx = std::min(endRotIdx, last);
  int rot=0, sidePrev=0;
  for(int i=orStart;i<=endRotIdx;++i){
    int side = (sc.Close[i]>vwap?1:(sc.Close[i]<vwap?-1:0));
    if(sidePrev!=0 && side!=0 && side!=sidePrev) rot++;
    if(side!=0) sidePrev=side;
  }

  // opening position vs prior day
  bool priorKnown = (Prior_VAL!=0 && Prior_VAH!=0);
  double openPx = sc.Open[p_SessionStartIdx];
  bool openInsideValue = (priorKnown && openPx>=Prior_VAL && openPx<=Prior_VAH);
  bool openOOR = false;
  if(priorKnown){
    double prevHi=sc.High[p_SessionStartIdx-1], prevLo=sc.Low[p_SessionStartIdx-1];
    openOOR = !(openPx>=prevLo && openPx<=prevHi);
  }

  int openType = 0; double conf=0.0;
  if (p_HadORBreak && heldMins >= 10){ openType=1; conf=std::min(1.0,0.5+heldMins/30.0); } // Drive
  else if (priorKnown && openInsideValue && rot>=2){ openType=4; conf=0.6 + std::min(0.2, rot*0.05); } // Auction In-Range
  else if (priorKnown && openOOR){ openType=5; conf=0.55; } // Auction Out-of-Range
  else if (p_HadORBreak && heldMins < 10){
    bool crossedOpp=false;
    if(p_ORFirstBreakIdx>0){
      bool breakUp = sc.High[p_ORFirstBreakIdx] > orHi;
      bool breakDn = sc.Low[p_ORFirstBreakIdx]  < orLo;
      for(int i=p_ORFirstBreakIdx;i<=last;++i){
        if(breakUp && sc.Low[i] < orLo){ crossedOpp=true; break; }
        if(breakDn && sc.High[i] > orHi){ crossedOpp=true; break; }
      }
    }
    if(crossedOpp){ openType=3; conf=0.6; } // Rejection-Reverse
    else          { openType=2; conf=0.6; } // Test-Drive
  }
  // hysteresis
  double dP = fabs(conf - p_LastOpenProb);
  if (dP < 0.15) conf = p_LastOpenProb; else p_LastOpenProb = (float)conf;
  if (openType != 0){
    if (p_OpenTypeSticky == 0){ p_OpenTypeSticky = openType; p_OpenTypeProb = (float)conf; }
    else if (openType != p_OpenTypeSticky && conf >= p_OpenTypeProb + 0.15f){ p_OpenTypeSticky=openType; p_OpenTypeProb=(float)conf; }
    else { p_OpenTypeProb = std::max(p_OpenTypeProb, (float)conf); }
  }

  // ----- WRITE FEEDS -----
  const std::string base = JoinPath(sc.DataFilesFolder(), In_OutDir.GetString());
  CreateDirectoryA(base.c_str(), NULL);

  // timestamps
  SCDateTime now = NowDT(sc);
  int date = sc.GetTradingDayDate(now);
  int time = sc.SecondsSinceStartOfDay(now); // seconds since midnight (keeps 'timestamp' numeric)

  // feed_bins.csv
  {
    std::ostringstream oss;
    oss << "timestamp,date,time,grid_ticks,tick_size,price_ticks,volume\n";
    for(size_t i=0;i<px.size();++i){
      oss << time << "," << date << "," << time << "," << grid << "," << sc.TickSize
          << "," << (long long)px[i] << "," << (long long)vv[i] << "\n";
    }
    AtomicWrite(base + "\\feed_bins.csv", oss.str());
  }

  // quick ADR20
  double ADR20=0.0; int days=0, idx=p_SessionStartIdx-1;
  while(idx>0 && days<20){
    int d = sc.GetTradingDayDate(sc.BaseDateTimeIn[idx]);
    int j=idx; double hi=sc.High[idx], lo=sc.Low[idx];
    while(j>0 && sc.GetTradingDayDate(sc.BaseDateTimeIn[j])==d){ hi=std::max(hi,sc.High[j]); lo=std::min(lo,sc.Low[j]); --j; }
    ADR20 += (hi-lo); days++; idx=j;
  }
  if(days>0) ADR20/=days;

  // feed_features.csv (now includes open/day probabilities)
  {
    double prevHi = (p_SessionStartIdx>0)? sc.High[p_SessionStartIdx-1]:0.0;
    double prevLo = (p_SessionStartIdx>0)? sc.Low[p_SessionStartIdx-1]:0.0;
    double gap = (ADR20>0? (openPx - sc.Close[p_SessionStartIdx-1]) / ADR20 : 0.0);
    (void)gap; // currently unused, keep for later

    // crude day-type proxy probs (keep simple & bounded)
    double p_day_trend   = std::min(0.95, std::max(0.05, 0.3*(heldMins/10.0) + (rot<=1?0.2:0.0) + (openOOR?0.2:0.0)));
    double p_day_balance = std::min(0.95, std::max(0.05, 0.5*(openInsideValue?1.0:0.0) + (rot>=2?0.3:0.0)));
    double norm = std::max(1e-6, p_day_trend + p_day_balance);
    p_day_trend/=norm; p_day_balance/=norm;

    std::ostringstream oss;
    oss << "timestamp,date,time,grid_ticks,tick_size,session_open,prior_poc,prior_val,prior_vah,prev_high,prev_low,adr20,"
           "or_low,or_high,sigma_bins,profile_poc_fallback,va_low_fallback,va_high_fallback,"
           "open_type_id,open_type_prob,day_trend_prob,day_balance_prob\n";
    oss << time << "," << date << "," << time << "," << grid << "," << sc.TickSize << ","
        << openPx << "," << Prior_POC << "," << Prior_VAL << "," << Prior_VAH << ","
        << prevHi << "," << prevLo << "," << ADR20 << ","
        << orLo << "," << orHi << "," << p_SigmaActive << ","
        << POC_fallback << "," << VAlo_fallback << "," << VAhi_fallback << ","
        << p_OpenTypeSticky << "," << p_OpenTypeProb << ","
        << p_day_trend << "," << p_day_balance << "\n";
    AtomicWrite(base + "\\feed_features.csv", oss.str());
  }

  // ----- READ RESULTS & DRAW -----
  double postPOC=NAN, postVAL50=NAN, postVAH50=NAN, lammax=NAN; double wbar[4]={0,0,0,0};
  {
    std::string path = base + "\\results.csv";
    std::ifstream fin(path.c_str());
    if(fin.good()){
      std::string line, lastline; std::getline(fin,line); // header
      while(std::getline(fin,line)) if(!line.empty()) lastline=line; fin.close();
      if(!lastline.empty()){
        std::stringstream ss(lastline); std::string tok; std::vector<std::string> c;
        while(std::getline(ss,tok,',')) c.push_back(tok);
        // header: timestamp,post_poc,post_va_low,post_va_high,post_poc_lo50,post_poc_hi50,post_va_low50,post_va_high50,lam_max,w1,w2,w3,w4
        if(c.size()>=13){
          postPOC   = atof(c[1].c_str());
          postVAL50 = atof(c[6].c_str()); // low50
          postVAH50 = atof(c[7].c_str()); // high50
          lammax    = atof(c[8].c_str());
          for(int k=0;k<4;++k) wbar[k] = atof(c[9+k].c_str());
          for(int i=sc.UpdateStartIndex;i<=last;++i){
            SG_PostPOC[i]=postPOC;
            SG_PostVA_Dn[i]=atof(c[2].c_str());
            SG_PostVA_Up[i]=atof(c[3].c_str());
            SG_LamMax[i]=lammax;
          }
          // credible band shading (50%)
          if(In_DrawBands.GetYesNo()){
            s_UseTool band; band.Clear(); band.ChartNumber=sc.ChartNumber; band.AddMethod=UTAM_ADD_OR_ADJUST;
            band.DrawingType=DRAWING_RECTANGLEHIGHLIGHT; band.Color=RGB(0,160,160); band.TransparencyLevel=88;
            band.Region=sc.GraphRegion; band.BeginIndex=p_SessionStartIdx; band.EndIndex=last;
            band.BeginValue=postVAL50; band.EndValue=postVAH50; sc.UseTool(band);
          }
        }
      }
    }
  }

  // ----- Flow markers at VAH/VAL/POC -----
  auto mark_flow = [&](double priceLevel, const char* labelBase){
    if(!(priceLevel>0)) return;
    int lookback = std::min(In_FlowLookback.GetInt(), last - p_SessionStartIdx);
    double bidVol=0, askVol=0;
    for(int bi=last - lookback; bi<=last; ++bi){
      int cnt = sc.VolumeAtPriceForBars->GetSizeAtBarIndex(bi);
      s_VolumeAtPriceV2* e=nullptr;
      for(int vi=0; vi<cnt; ++vi){
        sc.VolumeAtPriceForBars->GetVAPElementAtIndex(bi,vi,&e);
        if(!e) continue;
        double px = sc.TickSize * (double)e->PriceInTicks;
        if (fabs(px - priceLevel) <= sc.TickSize * 0.5){
          bidVol += (double)e->BidVolume;
          askVol += (double)e->AskVolume;
        }
      }
    }
    double thr = std::max(1.0, (double)In_ImbalanceThresh.GetFloat());
    const char* tag = nullptr;
    if(askVol/bidVol >= thr) tag = "SWEEP↑";
    else if(bidVol/askVol >= thr) tag = "SWEEP↓";
    else{
      // absorption: high combined vol but price not moving away much (last 5 bars range small)
      double hi=sc.High[last-5], lo=sc.Low[last-5];
      for(int i=last-5;i<=last;++i){ hi=std::max(hi,sc.High[i]); lo=std::min(lo,sc.Low[i]); }
      if((askVol+bidVol) > 0 && (hi-lo) <= 2*sc.TickSize) tag = (sc.LastTradePrice>=priceLevel? "ABSORB↑":"ABSORB↓");
    }
    if(tag){
      s_UseTool t; t.Clear(); t.ChartNumber=sc.ChartNumber; t.AddMethod=UTAM_ADD_OR_ADJUST;
      t.DrawingType=DRAWING_TEXT; t.Color=RGB(255,215,0); t.FontSize=10; t.Region=sc.GraphRegion;
      SCString txt; txt.Format("%s @ %s", tag, labelBase);
      t.Text = txt; t.BeginIndex=last; t.BeginValue=priceLevel; sc.UseTool(t);
    }
  };
  // Use posterior VA if available, else fallback
  double levPOC = (postPOC==postPOC? postPOC : POC_fallback);
  double levVAL = (postVAL50==postVAL50? postVAL50 : VAlo_fallback);
  double levVAH = (postVAH50==postVAH50? postVAH50 : VAhi_fallback);
  mark_flow(levVAL, "VAL"); mark_flow(levPOC, "POC"); mark_flow(levVAH, "VAH");

  // ----- HUD (posterior + weights bar) -----
  {
    auto bar = [](double x)->SCString{
      int n=10, f=(int)std::round(std::max(0.0,std::min(1.0,x))*n);
      SCString s("["); for(int i=0;i<n;++i) s += (i<f? "■":"□"); s += "]"; return s;
    };
    s_UseTool hud; hud.Clear(); hud.ChartNumber=sc.ChartNumber; hud.AddMethod=UTAM_ADD_OR_ADJUST;
    hud.DrawingType=DRAWING_TEXT; hud.Color=RGB(255,255,255); hud.FontSize=10; hud.Region=sc.GraphRegion;
    SCString text;
    text.Format("Posterior POC %.2f | VA50 [%.2f, %.2f] | λmax %.2f | w: S[%s] M[%s] L[%s]",
      levPOC, levVAL, levVAH, SG_LamMax[last],
      bar(wbar[0]).GetChars(), bar(wbar[1]).GetChars(), bar(wbar[2]).GetChars());
    hud.Text = text; hud.BeginIndex=last; hud.BeginValue=sc.LastTradePrice; sc.UseTool(hud);
  }
}
