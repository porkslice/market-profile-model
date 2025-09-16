# pythonside_marketprofilemodel.py — Priors from open/day, MH, bands, weights
import os, time, math
import numpy as np
import pandas as pd

try:
    import pywt
except ImportError:
    pywt = None

try:
    from sklearn.covariance import LedoitWolf
    HAVE_SK = True
except Exception:
    HAVE_SK = False

DATA_DIR = r"C:\SierraChart\Data\bayes"   # <-- set to your Sierra Data Files Folder\bayes

# ---------- IO ----------
def read_latest():
    bins = pd.read_csv(os.path.join(DATA_DIR,"feed_bins.csv"))
    feat = pd.read_csv(os.path.join(DATA_DIR,"feed_features.csv"))
    if bins.empty or feat.empty: return None, None
    ts = bins["timestamp"].max()
    bins = bins[bins["timestamp"]==ts].copy()
    return bins.sort_values("price_ticks"), feat.iloc[-1].to_dict()

def write_results(res):
    path = os.path.join(DATA_DIR,"results.csv")
    cols = ["timestamp","post_poc","post_va_low","post_va_high",
            "post_poc_lo50","post_poc_hi50","post_va_low50","post_va_high50",
            "lam_max","w1","w2","w3","w4"]
    df = pd.DataFrame([[res[k] for k in cols]], columns=cols)
    if not os.path.exists(path):
        df.to_csv(path, index=False)
    else:
        df.to_csv(path, mode="a", header=False, index=False)

# ---------- profile helpers ----------
def gaussian_smooth(x, sigma_bins):
    if sigma_bins<=0: return np.asarray(x,float)
    x=np.asarray(x,float); r=int(math.ceil(3*sigma_bins))
    k=np.arange(-r,r+1); w=np.exp(-(k*k)/(2*sigma_bins*sigma_bins)); w/=w.sum()
    return np.convolve(x,w,mode="same")

def poc_and_va(prob, prices, q=0.70):
    p=np.clip(np.asarray(prob,float),1e-20,None); p/=p.sum()
    idx=int(np.argmax(p)); L=R=idx; acc=p[idx]; n=p.size
    while acc<q and (L>0 or R+1<n):
        left=p[L-1] if L>0 else -1; right=p[R+1] if R+1<n else -1
        if right>=left and R+1<n: R+=1; acc+=p[R]
        elif L>0: L-=1; acc+=p[L]
        else: break
    return prices[idx], prices[L], prices[R], idx, L, R

def swt_components(curve, levels=3, wavelet="db4"):
    x=np.asarray(curve,float)
    if pywt is None:
        # Fallback: multiscale by smoothing differences
        sigs=[1,2,4][:levels]; comps=[]; prev=x.copy()
        for s in sigs:
            sm=gaussian_smooth(prev,s); comps.append(np.clip(prev-sm,0,None)); prev=sm
        comps.append(np.clip(prev,0,None)); comps=[c for c in comps if c.sum()>0]
    else:
        n=x.size; pow2=1<<(n-1).bit_length(); pad=pow2-n
        xp=np.pad(x,(0,pad),mode="edge")
        coeffs=pywt.swt(xp,wavelet,level=levels,trim_approx=True)
        comps=[np.abs(cD)[:n] for (cA,cD) in coeffs]  # details: small→mid
        comps.append(np.abs(coeffs[-1][0])[:n])       # approx: large
    # normalise each to sum 1
    out=[]
    for c in comps[:levels]:
        s=c.sum(); out.append(c/s if s>0 else np.ones_like(c)/len(c))
    return out[:levels]

# ---------- priors from context ----------
def map_priors(feat, S, base_sigma=2.0, sigma_bounds=(1.0,4.0)):
    # Open-type map
    ot_id  = int(feat.get("open_type_id",0))
    ot_prob= float(feat.get("open_type_prob",0.0))
    day_tr = float(feat.get("day_trend_prob",0.5))
    day_bal= float(feat.get("day_balance_prob",0.5))

    # weights prior (Dirichlet α) over S scales: [small, mid, large,...]
    alpha = np.ones(S)
    # Drive / Test-Drive bias small-scale
    drive_like = 1.0 if ot_id in (1,2) else 0.0
    auction_in = 1.0 if ot_id==4 else 0.0

    # blend open/day with caps (opening weight ≤ 0.4)
    w_open = min(0.4, max(0.0, ot_prob))
    w_day  = 1.0 - w_open
    trend_strength   = w_open*drive_like*ot_prob + w_day*day_tr
    balance_strength = w_open*auction_in*ot_prob + w_day*day_bal

    # allocate mass
    if S>=3:
        alpha[0] = 1.0 + 4.0*trend_strength       # small
        alpha[1] = 1.0 + 4.0*balance_strength     # mid
        alpha[2] = 1.0 + 2.0*(1.0 - abs(trend_strength-balance_strength))  # large stays flexible

    # sigma prior center
    sigma_center = base_sigma * (1.0 + 0.5*trend_strength - 0.3*balance_strength)
    sigma_center = float(np.clip(sigma_center, *sigma_bounds))

    priors = {"alpha": alpha, "sigma_center": sigma_center, "sigma_sd": 0.5}
    return priors

# ---------- MH ----------
rng = np.random.default_rng(123)

def dir_rw(w, kappa=220.0):
    a=np.maximum(1e-6, w*kappa)
    return rng.dirichlet(a)

def mh_fit(counts, prices, comps, sigma0=2.0, sigma_bounds=(1.0,4.0),
           iters=1200, burn=400, kappa=220.0, priors=None):
    K=len(prices); S=len(comps)
    comps=[np.clip(c,1e-20,None)/max(1e-20,c.sum()) for c in comps]
    counts=np.asarray(counts,float); N=counts.sum()
    w=np.ones(S)/S; sigma=float(np.clip(sigma0, *sigma_bounds))

    alpha=np.ones(S); sig_center=sigma; sig_sd=0.5
    if priors:
        alpha=np.asarray(priors.get("alpha",alpha),float)
        sig_center=float(priors.get("sigma_center",sig_center))
        sig_sd=float(priors.get("sigma_sd",sig_sd))

    def build_p(w,sig):
        mix=np.zeros(K)
        for s in range(S): mix += w[s]*comps[s]
        mix=gaussian_smooth(mix,sig); mix=np.clip(mix,1e-20,None); mix/=mix.sum()
        return mix

    def logpost(w,sig):
        p=build_p(w,sig)
        ll=float(np.sum(counts*np.log(p)))
        lpw=float(np.sum((alpha-1.0)*np.log(np.clip(w,1e-20,None))))
        lps=-0.5*((np.log(sig)-np.log(sig_center))/sig_sd)**2
        return ll+lpw+lps

    lp=logpost(w,sigma)
    W=[]; Sg=[]; accw=accs=0
    for it in range(iters):
        w2=dir_rw(w,kappa)
        lp2=logpost(w2,sigma)
        if np.log(rng.random()) < (lp2-lp): w=w2; lp=lp2; accw+=1
        s2=float(np.clip(np.exp(np.log(sigma)+rng.normal(0,0.15)), *sigma_bounds))
        lp3=logpost(w,s2)
        if np.log(rng.random()) < (lp3-lp): sigma=s2; lp=lp3; accs+=1
        if it>=burn: W.append(w.copy()); Sg.append(sigma)
    W=np.array(W); Sg=np.array(Sg)
    w_mean=W.mean(axis=0) if W.size else np.ones(S)/S
    sig_med=float(np.median(Sg)) if Sg.size else sigma
    # draws to get bands
    idxs=np.linspace(0, max(0,W.shape[0]-1), min(200,max(1,W.shape[0])) ,dtype=int)
    pocL=[]; vaL=[]; vaH=[]
    for j in idxs:
        p=build_p(W[j] if W.size else w, Sg[j] if Sg.size else sigma)
        poc, vl, vh, *_ = poc_and_va(p, prices, 0.70)
        pocL.append(poc); vaL.append(vl); vaH.append(vh)
    return {
        "w_mean": w_mean.tolist(),
        "sigma_median": sig_med,
        "poc_median": float(np.median(pocL)),
        "poc_lo50": float(np.percentile(pocL,25)),
        "poc_hi50": float(np.percentile(pocL,75)),
        "va_low_median": float(np.median(vaL)),
        "va_low50": float(np.percentile(vaL,25)),
        "va_high_median": float(np.median(vaH)),
        "va_high50": float(np.percentile(vaH,75)),
    }

# ---------- covariance ----------
class OnlineCov:
    def __init__(self, hl_sec=300):
        self.lam = math.exp(-math.log(2)/hl_sec)
        self.mu=None; self.S=None
    def update(self, x):
        lam=self.lam; x=np.asarray(x,float)
        if self.mu is None:
            self.mu=x.copy(); self.S=np.eye(x.size)*1e-6
        else:
            self.mu = lam*self.mu + (1-lam)*x
            dx = x - self.mu
            self.S  = lam*self.S + (1-lam)*np.outer(dx,dx)
        Ssh=self.S
        if HAVE_SK:
            try:
                Ssh = LedoitWolf().fit(self.S).covariance_
            except Exception: pass
        try:
            return float(np.linalg.eigvalsh(Ssh).max())
        except Exception:
            return float("nan")

ocov = OnlineCov(hl_sec=300)

# ---------- main loop ----------
def process_once():
    bins, feat = read_latest()
    if bins is None: return
    df = bins[["price_ticks","volume"]].groupby("price_ticks",as_index=False).sum().sort_values("price_ticks")
    prices = df["price_ticks"].values.astype(np.int64) * float(bins["tick_size"].iloc[0])
    counts = df["volume"].values.astype(float)
    mass = counts.sum()
    if mass <= 0: return
    curve = counts / mass

    # components
    S=3
    comps = swt_components(curve, levels=S, wavelet="db4")

    # priors from context
    priors = map_priors(feat, S=S, base_sigma=float(feat.get("sigma_bins",2.0)), sigma_bounds=(1.0,4.0))

    # MH
    post = mh_fit(counts, prices, comps,
                  sigma0=float(feat.get("sigma_bins",2.0)),
                  sigma_bounds=(1.0,4.0),
                  iters=1200, burn=400, kappa=220.0, priors=priors)

    # covariance input (weights + VA width + crude tails)
    va_w = post["va_high_median"] - post["va_low_median"]
    tail_up = max(0.0, prices.max() - post["va_high_median"])
    tail_dn = max(0.0, post["va_low_median"] - prices.min())
    x = np.array([*post["w_mean"], va_w, tail_up, tail_dn], float)
    lam_max = ocov.update(x)

    res = {
        "timestamp": feat["timestamp"],
        "post_poc": post["poc_median"],
        "post_va_low": post["va_low_median"],
        "post_va_high": post["va_high_median"],
        "post_poc_lo50": post["poc_lo50"],
        "post_poc_hi50": post["poc_hi50"],
        "post_va_low50": post["va_low50"],
        "post_va_high50": post["va_high50"],
        "lam_max": lam_max,
        "w1": (post["w_mean"] + [0.0])[0],
        "w2": (post["w_mean"] + [0.0])[1],
        "w3": (post["w_mean"] + [0.0])[2],
        "w4": (post["w_mean"] + [0.0,0.0])[3]  # pad
    }
    write_results(res)

def main():
    print("Sidecar watching", DATA_DIR)
    last_mtime = 0
    while True:
        try:
            m1 = os.path.getmtime(os.path.join(DATA_DIR,"feed_bins.csv")) if os.path.exists(os.path.join(DATA_DIR,"feed_bins.csv")) else 0
            m2 = os.path.getmtime(os.path.join(DATA_DIR,"feed_features.csv")) if os.path.exists(os.path.join(DATA_DIR,"feed_features.csv")) else 0
            sig = max(m1,m2)
            if sig > last_mtime:
                process_once()
                last_mtime = sig
        except KeyboardInterrupt:
            break
        except Exception as e:
            print("Error:", e)
        time.sleep(1.0)

if __name__ == "__main__":
    main()
