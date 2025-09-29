// marketprofile_bayes_inproc.cpp
#define NOMINMAX
#include "sierrachart.h"
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <limits>

#ifdef min
  #undef min
#endif
#ifdef max
  #undef max
#endif

SCDLLName("Bayesian Market-Profile – In-Process (v5.1)")

// ------------------ helpers ------------------
static constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

static inline SCDateTime AddMinutes(const SCDateTime& dt, int minutes){
  return dt + (minutes / 1440.0);
}

static inline int Round01To10Blocks(double x){
  if (x <= 0) return 0;
  if (x >= 1) return 10;
  return (int)std::floor(x * 10.0 + 0.5);
}

static std::vector<int> ParseIntList(const SCString& s){
  std::vector<int> out; out.reserve(4);
  std::string str = s.GetChars();
  std::stringstream ss(str);
  std::string tok;
  while(std::getline(ss, tok, ',')){
    if(!tok.empty()){
      try{ out.push_back(std::stoi(tok)); }catch(...){}
    }
  }
  if(out.empty()) out.push_back(2);
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

// smooth (gaussian over index)
static std::vector<double> SmoothBins(const std::vector<double>& x, double sigma_bins){
  if(x.empty()) return x;
  int n=(int)x.size();
  std::vector<double> y(n,0.0);
  int R=(int)std::ceil(3.0*sigma_bins);
  double denom = 2.0*sigma_bins*sigma_bins;
  for(int i=0;i<n;++i){
    double acc=0.0, wsum=0.0;
    int L = (i-R>0)?(i-R):0;
    int U = (i+R<n)?(i+R):(n-1);
    for(int j=L;j<=U;++j){
      double d=j-i;
      double w=std::exp(-(d*d)/denom);
      acc+=w*x[j]; wsum+=w;
    }
    y[i] = (wsum>0.0 ? acc/wsum : x[i]);
  }
  return y;
}

// normalize to prob-mass
static void NormalizeProb(std::vector<double>& p){
  double s=0.0;
  for(double v:p) if(v>0.0) s+=v;
  if(s<=0.0){ for(double&v:p) v=0.0; return; }
  for(double& v:p){ v = (v>0.0? v/s : 0.0); if(v<1e-20) v=1e-20; }
}

// POC + Value Area around POC covering `pct` of mass (0.5..0.99)
static void PocVaFromProb(const std::vector<int64_t>& px_ticks,
                          const std::vector<double>& prob,
                          double tickSize,
                          double pct,
                          double& pocPrice, double& vaLo, double& vaHi)
{
  pocPrice = vaLo = vaHi = 0.0;
  if(px_ticks.empty() || prob.empty() || px_ticks.size()!=prob.size()) return;

  double target = std::max(0.50, std::min(0.99, pct));
  size_t n = prob.size();
  size_t pocI = (size_t)(std::max_element(prob.begin(),prob.end()) - prob.begin());
  pocPrice = tickSize * (double)px_ticks[pocI];

  size_t L=pocI, R=pocI; double acc=prob[pocI];
  auto safe = [](double x){ return x>0.0? x:0.0; };
  while(acc<target && (L>0 || R+1<n)){
    double left  = (L>0)?prob[L-1]:-1.0;
    double right = (R+1<n)?prob[R+1]:-1.0;
    if(right>=left && R+1<n){ R++; acc+=safe(right); }
    else if(L>0){ L--; acc+=safe(left); }
    else break;
  }
  vaLo = tickSize * (double)px_ticks[L];
  vaHi = tickSize * (double)px_ticks[R];
}

// entropy-ish shape (λ ≈ exp(H0 - H))
static double EntropyLambda(const std::vector<double>& p){
  const double eps=1e-20;
  int n = (int)p.size();
  if(n<=1) return 1.0;
  double H=0.0;
  for(double v:p){
    double q = (v>eps? v:eps);
    H -= q*std::log(q);
  }
  double H0 = std::log((double)n);
  return std::exp(H0 - H); // 1 = uniform, grows as profile peaks
}

// quantile from a discrete pmf over price bins
static double PMFQuantile(const std::vector<int64_t>& px_ticks,
                          const std::vector<double>& prob,
                          double tickSize,
                          double q) // 0..1
{
  if(px_ticks.empty() || prob.empty() || px_ticks.size()!=prob.size()) return 0.0;
  q = std::max(0.0, std::min(1.0, q));
  double c=0.0;
  for(size_t i=0;i<prob.size();++i){
    c += prob[i];
    if(c >= q) return tickSize * (double)px_ticks[i];
  }
  return tickSize * (double)px_ticks.back();
}

// Bowley (quartile) skew: (Q3 + Q1 - 2*median) / (Q3 - Q1)
static double BowleySkew(const std::vector<int64_t>& px_ticks,
                         const std::vector<double>& prob,
                         double tickSize)
{
  if(px_ticks.size()<3 || prob.size()<3 || px_ticks.size()!=prob.size()) return 0.0;
  double q1 = PMFQuantile(px_ticks, prob, tickSize, 0.25);
  double q2 = PMFQuantile(px_ticks, prob, tickSize, 0.50);
  double q3 = PMFQuantile(px_ticks, prob, tickSize, 0.75);
  double denom = (q3 - q1);
  if(std::fabs(denom) < 1e-12) return 0.0;
  return (q3 + q1 - 2.0*q2) / denom;
}

// collect VAP for a grid
static void BuildGridVAP(SCStudyInterfaceRef sc, int startIdx, int lastIdx, int grid,
                         std::vector<int64_t>& px, std::vector<double>& vol)
{
  std::map<int64_t,double> vap;
  s_VolumeAtPriceV2* e=nullptr;
  for(int bi=startIdx; bi<=lastIdx; ++bi){
    int cnt = sc.VolumeAtPriceForBars->GetSizeAtBarIndex(bi);
    for(int vi=0; vi<cnt; ++vi){
      sc.VolumeAtPriceForBars->GetVAPElementAtIndex(bi,vi,&e);
      if(!e) continue;
      int64_t gbin = (e->PriceInTicks / grid) * (int64_t)grid;
      vap[gbin] += (double)e->Volume;
    }
  }
  px.clear(); vol.clear(); px.reserve(vap.size()); vol.reserve(vap.size());
  for(auto& kv : vap){ px.push_back(kv.first); vol.push_back(kv.second); }
}

// ------------------ study ------------------
SCSFExport scsf_BayesProfile_InProcess(SCStudyInterfaceRef sc)
{
  // Inputs
  SCInputRef In_GridsCsv        = sc.Input[0];
  SCInputRef In_SigmaBins       = sc.Input[1];
  SCInputRef In_SigmaMin        = sc.Input[2];
  SCInputRef In_SigmaMax        = sc.Input[3];
  SCInputRef In_VAPct           = sc.Input[4];
  SCInputRef In_OR_Min          = sc.Input[5];
  SCInputRef In_RefreshSec      = sc.Input[6];
  SCInputRef In_DrawPriorRefs   = sc.Input[7];
  SCInputRef In_FlowLookback    = sc.Input[8];
  SCInputRef In_ImbalanceThresh = sc.Input[9];
  SCInputRef In_DrawBands       = sc.Input[10];
  SCInputRef In_UICompact       = sc.Input[11];

  // Subgraphs
  SCSubgraphRef SG_POC       = sc.Subgraph[0];
  SCSubgraphRef SG_VA_Up     = sc.Subgraph[1];
  SCSubgraphRef SG_VA_Dn     = sc.Subgraph[2];
  SCSubgraphRef SG_PostPOC   = sc.Subgraph[3];
  SCSubgraphRef SG_PostVA_Up = sc.Subgraph[4];
  SCSubgraphRef SG_PostVA_Dn = sc.Subgraph[5];

  // Persistent state
  int&   p_SessionStartIdx  = sc.GetPersistentInt(1);
  int&   p_OR_EndIdx        = sc.GetPersistentInt(2);
  int&   p_LastRefreshEpoch = sc.GetPersistentInt(3);
  int&   p_DayKey           = sc.GetPersistentInt(4);
  float& p_SigmaActive      = sc.GetPersistentFloat(1);

  double& Prior_VAH = sc.GetPersistentDouble(1);
  double& Prior_VAL = sc.GetPersistentDouble(2);
  double& Prior_POC = sc.GetPersistentDouble(3);

  int&   p_OpenTypeSticky   = sc.GetPersistentInt(5);
  float& p_OpenTypeProb     = sc.GetPersistentFloat(2);
  float& p_LastOpenProb     = sc.GetPersistentFloat(3);
  int&   p_HadORBreak       = sc.GetPersistentInt(6);
  int&   p_ORFirstBreakIdx  = sc.GetPersistentInt(7);

  int&   d_OR_Box           = sc.GetPersistentInt(20);
  int&   d_PriorVAL         = sc.GetPersistentInt(21);
  int&   d_PriorVAH         = sc.GetPersistentInt(22);
  int&   d_PriorPOC         = sc.GetPersistentInt(23);
  int&   d_CredBand         = sc.GetPersistentInt(24);
  int&   d_HUD              = sc.GetPersistentInt(25);
  int&   d_FlowVAL          = sc.GetPersistentInt(26);
  int&   d_FlowPOC          = sc.GetPersistentInt(27);
  int&   d_FlowVAH          = sc.GetPersistentInt(28);

  double& p_LastLam         = sc.GetPersistentDouble(10);
  double& p_LastBskew       = sc.GetPersistentDouble(11);

  if (sc.SetDefaults)
  {
    sc.GraphName = "Bayesian Market-Profile – In-Process (v5.1)";
    sc.AutoLoop = 0; sc.GraphRegion = 0;
    sc.MaintainVolumeAtPriceData = 1; sc.UpdateAlways = 1;

    In_GridsCsv.Name  = "Grids (CSV)";                 In_GridsCsv.SetString("1,2,4");
    In_SigmaBins.Name = "Smoothing Sigma (bins)";      In_SigmaBins.SetFloat(2.0f);
    In_SigmaMin.Name  = "Sigma Min (bins)";            In_SigmaMin.SetFloat(1.0f);
    In_SigmaMax.Name  = "Sigma Max (bins)";            In_SigmaMax.SetFloat(4.0f);
    In_VAPct.Name     = "Value Area Percent (0..1)";   In_VAPct.SetFloat(0.70f);
    In_OR_Min.Name    = "Opening Range Minutes";       In_OR_Min.SetInt(5); In_OR_Min.SetIntLimits(1,30);
    In_RefreshSec.Name= "Refresh Seconds";             In_RefreshSec.SetInt(8); In_RefreshSec.SetIntLimits(1,30);
    In_DrawPriorRefs.Name = "Draw Prior VAH/VAL/POC";  In_DrawPriorRefs.SetYesNo(true);
    In_FlowLookback.Name = "Flow Lookback Bars";       In_FlowLookback.SetInt(10); In_FlowLookback.SetIntLimits(3,50);
    In_ImbalanceThresh.Name = "Flow Imbalance Ratio";  In_ImbalanceThresh.SetFloat(1.8f);
    In_DrawBands.Name = "Draw VA 50% Band";            In_DrawBands.SetYesNo(true);
    In_UICompact.Name = "Compact HUD";                 In_UICompact.SetYesNo(false);

    SG_POC.Name="POC (fallback)"; SG_POC.DrawStyle=DRAWSTYLE_LINE; SG_POC.PrimaryColor=RGB(200,200,200); SG_POC.LineWidth=2;
    SG_VA_Up.Name="VA Upper (fallback)"; SG_VA_Up.DrawStyle=DRAWSTYLE_LINE; SG_VA_Up.PrimaryColor=RGB(0,200,200); SG_VA_Up.LineStyle=LINESTYLE_DASH;
    SG_VA_Dn.Name="VA Lower (fallback)"; SG_VA_Dn.DrawStyle=DRAWSTYLE_LINE; SG_VA_Dn.PrimaryColor=RGB(0,200,200); SG_VA_Dn.LineStyle=LINESTYLE_DASH;

    SG_PostPOC.Name="POC (posterior)"; SG_PostPOC.DrawStyle=DRAWSTYLE_LINE; SG_PostPOC.PrimaryColor=RGB(255,255,255); SG_PostPOC.LineWidth=2;
    SG_PostVA_Up.Name="VA Upper (posterior)"; SG_PostVA_Up.DrawStyle=DRAWSTYLE_LINE; SG_PostVA_Up.PrimaryColor=RGB(255,255,255); SG_PostVA_Up.LineStyle=LINESTYLE_DOT;
    SG_PostVA_Dn.Name="VA Lower (posterior)"; SG_PostVA_Dn.DrawStyle=DRAWSTYLE_LINE; SG_PostVA_Dn.PrimaryColor=RGB(255,255,255); SG_PostVA_Dn.LineStyle=LINESTYLE_DOT;

    return;
  }

  const int last = sc.ArraySize - 1; if (last < 10) return;

  // ---- New trading day detection + prior refs (unchanged) ----
  const int todayKey = sc.GetTradingDayDate(sc.BaseDateTimeIn[last]);
  bool is_new_td = (p_DayKey != todayKey);
  if (!is_new_td && last>0){
    int prevKey = sc.GetTradingDayDate(sc.BaseDateTimeIn[last-1]);
    is_new_td = (prevKey != todayKey);
  }
  if (is_new_td){
    p_DayKey = todayKey;
    p_SessionStartIdx = last;
    for (int i=last;i>=0;--i){
      if (sc.GetTradingDayDate(sc.BaseDateTimeIn[i]) != todayKey) break;
      p_SessionStartIdx = i;
    }
    SCDateTime orEnd = AddMinutes(sc.BaseDateTimeIn[p_SessionStartIdx], In_OR_Min.GetInt());
    p_OR_EndIdx = sc.GetContainingIndexForSCDateTime(sc.ChartNumber, orEnd);
    p_SigmaActive = (float)In_SigmaBins.GetFloat();

    p_OpenTypeSticky=0; p_OpenTypeProb=0.0f; p_LastOpenProb=0.0f;
    p_HadORBreak=0; p_ORFirstBreakIdx=-1;

    for (int i=p_SessionStartIdx; i<=last; ++i){
      SG_PostPOC[i]=kNaN; SG_PostVA_Dn[i]=kNaN; SG_PostVA_Up[i]=kNaN;
    }
    p_LastLam = 0.0; p_LastBskew = 0.0;

    Prior_VAH = Prior_VAL = Prior_POC = 0.0;
    if (In_DrawPriorRefs.GetYesNo()){
      int prevEnd=-1, prevStart=-1;
      for(int i=p_SessionStartIdx-1; i>=0; --i){
        if(prevEnd==-1) prevEnd=i;
        if(sc.GetTradingDayDate(sc.BaseDateTimeIn[i]) != todayKey){ prevStart=i+1; break; }
      }
      if(prevStart>=0 && prevEnd>prevStart){
        std::map<int64_t,double> vap; s_VolumeAtPriceV2* e=nullptr;
        for(int bi=prevStart; bi<=prevEnd; ++bi){
          int cnt = sc.VolumeAtPriceForBars->GetSizeAtBarIndex(bi);
          for(int vi=0; vi<cnt; ++vi){
            sc.VolumeAtPriceForBars->GetVAPElementAtIndex(bi,vi,&e);
            if(!e) continue; vap[e->PriceInTicks] += (double)e->Volume;
          }
        }
        if(!vap.empty()){
          std::vector<int64_t> px; std::vector<double> vv;
          double total=0; px.reserve(vap.size()); vv.reserve(vap.size());
          for(auto& kv:vap){ px.push_back(kv.first); vv.push_back(kv.second); total+=kv.second; }
          size_t pocI = (size_t)(std::max_element(vv.begin(),vv.end())-vv.begin());
          Prior_POC = sc.TickSize * (double)px[pocI];
          double target = total * In_VAPct.GetFloat(); size_t L=pocI, R=pocI; double acc=vv[pocI];
          while(acc<target && (L>0 || R+1<vv.size())){
            double left =(L>0)?vv[L-1]:-1.0; double right=(R+1<vv.size())?vv[R+1]:-1.0;
            if(right>=left && R+1<vv.size()){ R++; acc+=right; }
            else if(L>0){ L--; acc+=left; } else break;
          }
          Prior_VAL = sc.TickSize*(double)px[L];
          Prior_VAH = sc.TickSize*(double)px[R];
          s_UseTool t; t.Clear(); t.ChartNumber=sc.ChartNumber; t.AddMethod=UTAM_ADD_OR_ADJUST; t.Region=sc.GraphRegion;
          t.DrawingType=DRAWING_HORIZONTALLINE; t.LineStyle=LINESTYLE_DOT; t.LineWidth=1;
          t.Color=RGB(255,0,200); t.BeginValue=Prior_VAH; sc.UseTool(t); d_PriorVAH=t.LineNumber;
          t.BeginValue=Prior_VAL; sc.UseTool(t); d_PriorVAL=t.LineNumber;
          t.Color=RGB(200,200,200); t.LineStyle=LINESTYLE_DASH; t.BeginValue=Prior_POC; sc.UseTool(t); d_PriorPOC=t.LineNumber;
        }
      }
    }
  }

  // ---- Throttle ----
  const int refresh = In_RefreshSec.GetInt();
  const int nowSec = sc.BaseDateTimeIn[last].GetTimeInSeconds();
  if (nowSec - p_LastRefreshEpoch < refresh && sc.UpdateStartIndex!=0) return;
  p_LastRefreshEpoch = nowSec;

  // ---- OR box ----
  int orStart=p_SessionStartIdx; int orEnd=p_OR_EndIdx; if (orEnd > last) orEnd = last;
  double orHi=0, orLo=0;
  if(orEnd>orStart){
    orHi=sc.High[orStart]; orLo=sc.Low[orStart];
    for (int i = orStart; i <= orEnd; ++i){
      if (sc.High[i] > orHi) orHi = sc.High[i];
      if (sc.Low[i]  < orLo) orLo = sc.Low[i];
    }
    s_UseTool t; t.Clear(); t.ChartNumber=sc.ChartNumber; t.Region=sc.GraphRegion; t.AddMethod=UTAM_ADD_OR_ADJUST;
    t.DrawingType=DRAWING_RECTANGLEHIGHLIGHT; t.Color=RGB(70,130,180); t.TransparencyLevel=85;
    t.BeginIndex=orStart; t.EndIndex=orEnd; t.BeginValue=orLo; t.EndValue=orHi; t.LineNumber=d_OR_Box; sc.UseTool(t); d_OR_Box=t.LineNumber;
  }

  // ---- Opening type heuristics (unchanged) ----
  if (!p_HadORBreak && last > orEnd){
    for(int i=orEnd+1;i<=last;++i){
      if(sc.High[i] > orHi || sc.Low[i] < orLo){ p_HadORBreak=1; p_ORFirstBreakIdx=i; break; }
    }
  }
  double heldMins=0.0;
  if (p_HadORBreak && p_ORFirstBreakIdx>0){
    bool outsideNow = (sc.LastTradePrice > orHi || sc.LastTradePrice < orLo);
    if(outsideNow){
      int t1 = sc.BaseDateTimeIn[p_ORFirstBreakIdx].GetTimeInSeconds();
      int t2 = sc.BaseDateTimeIn[last].GetTimeInSeconds();
      heldMins = (t2 - t1) / 60.0;
    }
  }
  double vwap=0, pv=0, vol=0; for(int i=p_SessionStartIdx;i<=last;++i){ pv += sc.Volume[i]*sc.HLCAvg[i]; vol += sc.Volume[i]; }
  if(vol>0) vwap=pv/vol;
  SCDateTime halfHour = AddMinutes(sc.BaseDateTimeIn[orStart], 30);
  int endRotIdx = sc.GetContainingIndexForSCDateTime(sc.ChartNumber, halfHour);
  if (endRotIdx > last) endRotIdx = last;
  int rot=0, sidePrev=0;
  for(int i=orStart;i<=endRotIdx;++i){
    int side = (sc.Close[i]>vwap?1:(sc.Close[i]<vwap?-1:0));
    if(sidePrev!=0 && side!=0 && side!=sidePrev) rot++;
    if(side!=0) sidePrev=side;
  }
  bool priorKnown = (Prior_VAL!=0 && Prior_VAH!=0);
  double openPx = sc.Open[p_SessionStartIdx];
  bool openInsideValue = (priorKnown && openPx>=Prior_VAL && openPx<=Prior_VAH);
  bool openOOR=false; if(priorKnown){ double prevHi=sc.High[p_SessionStartIdx-1], prevLo=sc.Low[p_SessionStartIdx-1]; openOOR = !(openPx>=prevLo && openPx<=prevHi); }
  int openType = 0; double conf=0.0;
  if (p_HadORBreak && heldMins >= 10){ openType=1; conf=(0.5+heldMins/30.0); if (conf>1.0) conf=1.0; }
  else if (priorKnown && openInsideValue && rot>=2){ openType=4; conf=0.6 + (rot*0.05); if (conf>0.8) conf=0.8; }
  else if (priorKnown && openOOR){ openType=5; conf=0.55; }
  else if (p_HadORBreak && heldMins < 10){
    bool crossedOpp=false;
    if(p_ORFirstBreakIdx>0){
      bool breakUp = sc.High[p_ORFirstBreakIdx] > orHi;
      bool breakDn = sc.Low[p_ORFirstBreakIdx]  < orLo;
      for(int i=p_ORFirstBreakIdx;i<=last;++i){
        if(breakUp && sc.Low[i]  < orLo){ crossedOpp=true; break; }
        if(breakDn && sc.High[i] > orHi){ crossedOpp=true; break; }
      }
    }
    if(crossedOpp){ openType=3; conf=0.6; } else { openType=2; conf=0.6; }
  }
  double dP = std::fabs(conf - p_LastOpenProb);
  if (dP < 0.15) conf = p_LastOpenProb; else p_LastOpenProb = (float)conf;
  if (openType != 0){
    if (p_OpenTypeSticky == 0){ p_OpenTypeSticky = openType; p_OpenTypeProb = (float)conf; }
    else if (openType != p_OpenTypeSticky && conf >= p_OpenTypeProb + 0.15f){ p_OpenTypeSticky=openType; p_OpenTypeProb=(float)conf; }
    else { if ((float)conf > p_OpenTypeProb) p_OpenTypeProb = (float)conf; }
  }

  // ---------- Multi-grid build ----------
  // ADR20 → volatility ratio
  double ADR20=0.0; int days=0, idx=p_SessionStartIdx-1;
  while(idx>0 && days<20){
    int d = sc.GetTradingDayDate(sc.BaseDateTimeIn[idx]);
    int j=idx; double hi=sc.High[idx], lo=sc.Low[idx];
    while (j > 0 && sc.GetTradingDayDate(sc.BaseDateTimeIn[j]) == d){
      if (sc.High[j] > hi) hi = sc.High[j];
      if (sc.Low[j]  < lo) lo = sc.Low[j];
      --j;
    }
    ADR20 += (hi-lo); days++; idx=j;
  }
  if(days>0) ADR20/=days;

  double sessHi=sc.High[p_SessionStartIdx], sessLo=sc.Low[p_SessionStartIdx];
  for(int i=p_SessionStartIdx;i<=last;++i){ if(sc.High[i]>sessHi) sessHi=sc.High[i]; if(sc.Low[i]<sessLo) sessLo=sc.Low[i]; }
  double sessRange = std::max(1e-9, sessHi - sessLo);
  double volRatio = (ADR20>0.0? sessRange/ADR20 : 1.0);

  // sigma tracking
  float sMin=(float)In_SigmaMin.GetFloat();
  float sMax=(float)In_SigmaMax.GetFloat();
  float sNew=(float)In_SigmaBins.GetFloat();
  if (sNew < sMin) sNew = sMin; if (sNew > sMax) sNew = sMax;
  if (p_SigmaActive<=0) p_SigmaActive = sNew;
  else {
    float denom = (p_SigmaActive!=0.0f)? p_SigmaActive : 1e-6f;
    float rel = std::fabs(sNew - p_SigmaActive)/denom;
    if (rel >= 0.10f) p_SigmaActive = sNew;
  }

  std::vector<int> grids = ParseIntList(In_GridsCsv.GetString());
  if(grids.empty()) grids.push_back(2);

  struct GRes{ int g; double poc, vaL, vaH, vaL50, vaH50, lam, bskew; };
  std::vector<GRes> gres; gres.reserve(grids.size());

  for(int g : grids){
    std::vector<int64_t> px; std::vector<double> vv;
    BuildGridVAP(sc, p_SessionStartIdx, last, g, px, vv);
    if(px.size()<5) continue;

    std::vector<double> sm = SmoothBins(vv, p_SigmaActive);
    NormalizeProb(sm);

    double poc, vaL, vaH; PocVaFromProb(px, sm, sc.TickSize, In_VAPct.GetFloat(), poc, vaL, vaH);
    double L50,H50; PocVaFromProb(px, sm, sc.TickSize, 0.50, poc/*unused*/, L50, H50);
    double lam = EntropyLambda(sm);
    double bsk = BowleySkew(px, sm, sc.TickSize);

    gres.push_back({g,poc,vaL,vaH,L50,H50,lam,bsk});
  }

  double POC_fb=0, VAL_fb=0, VAH_fb=0;
  if(!gres.empty()){
    auto base = *std::min_element(gres.begin(), gres.end(), [](const GRes& a, const GRes& b){ return a.g < b.g; });
    POC_fb = base.poc; VAL_fb=base.vaL; VAH_fb=base.vaH;
  }
  for(int i=sc.UpdateStartIndex;i<=last;++i){
    SG_POC[i]=POC_fb; SG_VA_Dn[i]=VAL_fb; SG_VA_Up[i]=VAH_fb;
  }

  // weights S/M/L by volatility
  double wS=0.33,wM=0.34,wL=0.33;
  {
    double v = std::max(0.2, std::min(1.8, volRatio));
    auto lerp3 = [](double a1,double b1,double c1, double a2,double b2,double c2, double t){
      return std::tuple<double,double,double>(a1 + (a2-a1)*t, b1 + (b2-b1)*t, c1 + (c2-c1)*t);
    };
    if(v <= 0.75){
      double t = (v-0.5)/0.25; if(t<0) t=0; if(t>1) t=1;
      std::tie(wS,wM,wL) = lerp3(0.60,0.30,0.10, 0.40,0.45,0.15, t);
    } else if(v <= 1.25){
      double t = (v-1.0)/0.25; if(t<0) t=0; if(t>1) t=1;
      std::tie(wS,wM,wL) = lerp3(0.40,0.45,0.15, 0.25,0.50,0.25, t);
    } else {
      double t = (v-1.5)/0.3; if(t<0) t=0; if(t>1) t=1;
      std::tie(wS,wM,wL) = lerp3(0.25,0.50,0.25, 0.10,0.30,0.60, t);
    }
    double s=wS+wM+wL; if(s<=1e-9){ wS=wM=wL=1.0/3.0; } else { wS/=s; wM/=s; wL/=s; }
    std::sort(gres.begin(), gres.end(), [](const GRes&a,const GRes&b){ return a.g<b.g; });
  }

  double postPOC=0, postVAL=0, postVAH=0, postVAL50=0, postVAH50=0, lamMix=0, bskewMix=0;
  if(gres.size()==1){
    postPOC=gres[0].poc; postVAL=gres[0].vaL; postVAH=gres[0].vaH; postVAL50=gres[0].vaL50; postVAH50=gres[0].vaH50; lamMix=gres[0].lam; bskewMix=gres[0].bskew;
  } else if(gres.size()==2){
    double wA=wS+wM, wB=wL;
    postPOC = wA* (0.5*(gres[0].poc+gres[1].poc)) + wB*gres.back().poc;
    postVAL = wA* (0.5*(gres[0].vaL+gres[1].vaL)) + wB*gres.back().vaL;
    postVAH = wA* (0.5*(gres[0].vaH+gres[1].vaH)) + wB*gres.back().vaH;
    postVAL50 = wA* (0.5*(gres[0].vaL50+gres[1].vaL50)) + wB*gres.back().vaL50;
    postVAH50 = wA* (0.5*(gres[0].vaH50+gres[1].vaH50)) + wB*gres.back().vaH50;
    lamMix = wA*(0.5*(gres[0].lam+gres[1].lam)) + wB*gres.back().lam;
    bskewMix = wA*(0.5*(gres[0].bskew+gres[1].bskew)) + wB*gres.back().bskew;
  } else if(gres.size()>=3){
    postPOC  = wS*gres[0].poc   + wM*gres[1].poc   + wL*gres.back().poc;
    postVAL  = wS*gres[0].vaL   + wM*gres[1].vaL   + wL*gres.back().vaL;
    postVAH  = wS*gres[0].vaH   + wM*gres[1].vaH   + wL*gres.back().vaH;
    postVAL50= wS*gres[0].vaL50 + wM*gres[1].vaL50 + wL*gres.back().vaL50;
    postVAH50= wS*gres[0].vaH50 + wM*gres[1].vaH50 + wL*gres.back().vaH50;
    lamMix   = wS*gres[0].lam   + wM*gres[1].lam   + wL*gres.back().lam;
    bskewMix = wS*gres[0].bskew + wM*gres[1].bskew + wL*gres.back().bskew;
  }

  if(postPOC>0 && postVAL>0 && postVAH>0){
    for(int i=sc.UpdateStartIndex;i<=last;++i){
      SG_PostPOC[i]=postPOC; SG_PostVA_Dn[i]=postVAL; SG_PostVA_Up[i]=postVAH;
    }
    if(In_DrawBands.GetYesNo() && postVAL50>0 && postVAH50>0){
      s_UseTool band; band.Clear(); band.ChartNumber=sc.ChartNumber; band.AddMethod=UTAM_ADD_OR_ADJUST; band.Region=sc.GraphRegion;
      band.DrawingType=DRAWING_RECTANGLEHIGHLIGHT; band.Color=RGB(0,160,160); band.TransparencyLevel=88;
      band.BeginIndex=p_SessionStartIdx; band.EndIndex=last; band.BeginValue=postVAL50; band.EndValue=postVAH50; band.LineNumber=d_CredBand;
      sc.UseTool(band); d_CredBand=band.LineNumber;
    }
    p_LastLam   = lamMix;
    p_LastBskew = bskewMix;
  }

  // ---- Flow markers ----
  auto mark_flow = [&](double priceLevel, const char* labelBase, int& id_holder){
    if(!(priceLevel>0)) return;
    int lookback = In_FlowLookback.GetInt();
    int maxlb = last - p_SessionStartIdx; if (lookback > maxlb) lookback = maxlb; if (lookback < 0) lookback = 0;
    double bidVol=0, askVol=0; s_VolumeAtPriceV2* e=nullptr;
    for(int bi=last - lookback; bi<=last; ++bi){
      int cnt = sc.VolumeAtPriceForBars->GetSizeAtBarIndex(bi);
      for(int vi=0; vi<cnt; ++vi){
        sc.VolumeAtPriceForBars->GetVAPElementAtIndex(bi,vi,&e);
        if(!e) continue;
        double px = sc.TickSize * (double)e->PriceInTicks;
        if (std::fabs(px - priceLevel) <= sc.TickSize * 0.5){
          bidVol += (double)e->BidVolume; askVol += (double)e->AskVolume;
        }
      }
    }
    double thr = std::max(1.0, (double)In_ImbalanceThresh.GetFloat());
    const double eps = 1e-12; const char* tag = nullptr;
    if (askVol / (bidVol + eps) >= thr)      tag = "SWEEP\xE2\x86\x91";
    else if (bidVol / (askVol + eps) >= thr) tag = "SWEEP\xE2\x86\x93";
    else{
      int start = std::max(p_SessionStartIdx, last-5);
      double hi=sc.High[start], lo=sc.Low[start];
      for (int i = start; i <= last; ++i){ if (sc.High[i] > hi) hi=sc.High[i]; if (sc.Low[i] < lo) lo=sc.Low[i]; }
      if((askVol+bidVol) > 0 && (hi-lo) <= 2*sc.TickSize) tag = (sc.LastTradePrice>=priceLevel? "ABSORB\xE2\x86\x91":"ABSORB\xE2\x86\x93");
    }
    if(tag){
      s_UseTool t; t.Clear(); t.ChartNumber=sc.ChartNumber; t.AddMethod=UTAM_ADD_OR_ADJUST;
      t.Region=sc.GraphRegion; t.DrawingType=DRAWING_TEXT; t.Color=RGB(255,215,0); t.FontSize=10;
      SCString txt; txt.Format("%s @ %s", tag, labelBase); t.Text = txt; t.BeginIndex=last; t.BeginValue=priceLevel; t.LineNumber=id_holder;
      sc.UseTool(t); id_holder=t.LineNumber;
    }
  };
  const double levPOC = (SG_PostPOC[last]>0? SG_PostPOC[last] : POC_fb);
  const double levVAL = (SG_PostVA_Dn[last]>0? SG_PostVA_Dn[last] : VAL_fb);
  const double levVAH = (SG_PostVA_Up[last]>0? SG_PostVA_Up[last] : VAH_fb);
  mark_flow(levVAL,"VAL",d_FlowVAL); mark_flow(levPOC,"POC",d_FlowPOC); mark_flow(levVAH,"VAH",d_FlowVAH);

  // ---- HUD (now includes Bowley skew) ----
  {
    auto bar10 = [](double x)->SCString{ int b=Round01To10Blocks(x); SCString s("["); for(int i=0;i<10;++i) s += (i<b? "\xE2\x96\xA0":"\xE2\x96\xA1"); s += "]"; return s; };
    double shownS=0, shownM=0, shownL=0;
    if(gres.size()==1){ shownS=1; }
    else if(gres.size()==2){ shownS=wS+wM; shownL=wL; }
    else { shownS=wS; shownM=wM; shownL=wL; }

    const char* btag = (std::fabs(p_LastBskew) < 0.2 ? "Neutral" : (p_LastBskew > 0 ? "Right" : "Left"));

    SCString text;
    if (In_UICompact.GetYesNo()){
      text.Format("POC %.2f  VA[%.2f,%.2f]  VA50[%.2f,%.2f]  \xCE\xBB %.2f  Bskew %.2f (%s)  w S/M/L %.0f/%.0f/%.0f%%  Open:%d(%.0f%%)",
                  levPOC, levVAL, levVAH, postVAL50, postVAH50, p_LastLam, p_LastBskew, btag,
                  shownS*100.0, shownM*100.0, shownL*100.0, p_OpenTypeSticky, p_OpenTypeProb*100.0);
    } else {
      SCString wSbar = bar10(shownS), wMbar = bar10(shownM), wLbar = bar10(shownL);
      text.Format("Posterior POC %.2f\nVA [%.2f, %.2f]\nVA50 [%.2f, %.2f]\n\xCE\xBB %.2f\nBowley Skew %.2f (%s)\nWeights S:%s  M:%s  L:%s\nOpenType %d  P:%.0f%%",
                  levPOC, levVAL, levVAH, postVAL50, postVAH50, p_LastLam, p_LastBskew, btag,
                  wSbar.GetChars(), wMbar.GetChars(), wLbar.GetChars(), p_OpenTypeSticky, p_OpenTypeProb*100.0);
    }
    s_UseTool hud; hud.Clear(); hud.ChartNumber=sc.ChartNumber; hud.AddMethod=UTAM_ADD_OR_ADJUST;
    hud.Region=sc.GraphRegion; hud.DrawingType=DRAWING_TEXT; hud.Color=RGB(255,255,255); hud.FontSize=10;
    hud.Text=text; hud.BeginIndex=last; hud.BeginValue=sc.LastTradePrice; hud.LineNumber=d_HUD; sc.UseTool(hud); d_HUD=hud.LineNumber;
  }
}
