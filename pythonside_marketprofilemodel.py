# sidecar_main.py  — MH + SWT + covariance for Sierra bridge
import os, time, json, math, threading
import numpy as np
import pandas as pd

# Optional deps
try:
    import pywt
except ImportError:
    pywt = None

try:
    from sklearn.covariance import LedoitWolf
    _HAVE_SK = True
except Exception:
    _HAVE_SK = False

DATA_DIR = r"C:\SierraChart\Data\bayes"   # <-- CHANGE to your Sierra Data Files Folder\bayes

# ----------------------- utils -----------------------
def read_latest_bins_features():
    bins_path = os.path.join(DATA_DIR, "feed_bins.csv")
    feat_path = os.path.join(DATA_DIR, "feed_features.csv")
    if not (os.path.exists(bins_path) and os.path.exists(feat_path)):
        return None, None
    bins = pd.read_csv(bins_path)
    feats = pd.read_csv(feat_path)
    if bins.empty or feats.empty:
        return None, None
    # Use the single latest snapshot (largest timestamp)
    ts = bins["timestamp"].max()
    bins = bins[bins["timestamp"] == ts].copy()
    feat = feats.iloc[-1].to_dict()
    return bins, feat

def write_results(res):
    out_path = os.path.join(DATA_DIR, "results.csv")
    header = ["timestamp","post_poc","post_va_low","post_va_high",
              "post_poc_lo50","post_poc_hi50","post_va_low50","post_va_high50",
              "lam_max","w1","w2","w3","w4"]
    df = pd.DataFrame([[
        res["timestamp"], res["post_poc"], res["post_va_low"], res["post_va_high"],
        res["post_poc_lo50"], res["post_poc_hi50"], res["post_va_low50"], res["post_va_high50"],
        res["lam_max"], *res["weights_extended"]
    ]], columns=header)
    # append or create
    if not os.path.exists(out_path):
        df.to_csv(out_path, index=False)
    else:
        df.to_csv(out_path, mode="a", header=False, index=False)

# ------------------- profile helpers -------------------
def swt_components(curve, wavelet="db4", levels=3):
    """Return nonnegative per-scale components that sum to ~curve mass."""
    if pywt is None:
        # fallback: crude multi-scale via Gaussian smoothing differences
        comps = []
        sigmas = [1, 2, 4][:levels]
        prev = np.array(curve, float)
        for s in sigmas:
            sm = gaussian_smooth(prev, s)
            comp = np.clip(prev - sm, 0, None)
            comps.append(comp)
            prev = sm
        comps.append(np.clip(prev, 0, None))
        return comps[:levels]
    x = np.asarray(curve, float)
    # SWT needs length as power-of-two padding; pad reflect
    n = x.size
    pow2 = 1 << (n-1).bit_length()
    pad = pow2 - n
    xpad = np.pad(x, (0, pad), mode="edge")
    coeffs = pywt.swt(xpad, wavelet, level=levels, start_level=0, trim_approx=True)
    # coeffs: list of (cA_k, cD_k) pairs from level=1..levels
    comps = []
    recon = np.zeros_like(xpad)
    # Make each component nonnegative by using squared (energy-like) or abs
    for k, (cA, cD) in enumerate(coeffs, 1):
        comp = np.abs(cD)
        comps.append(comp[:n])
    # final approx as larger-scale
    comps.append(np.abs(coeffs[-1][0])[:n])
    # Normalize each component to sum 1 (when possible)
    for i in range(len(comps)):
        s = comps[i].sum()
        comps[i] = comps[i] / s if s > 0 else comps[i]
    return comps[:levels]

def gaussian_smooth(x, sigma_bins):
    if sigma_bins <= 0: return np.array(x, float)
    x = np.asarray(x, float)
    n = x.size
    r = int(math.ceil(3*sigma_bins))
    kidx = np.arange(-r, r+1)
    w = np.exp(-(kidx**2)/(2*sigma_bins**2))
    w /= w.sum()
    y = np.convolve(x, w, mode="same")
    return y

def poc_and_va(prob, prices, q=0.70):
    """POC index and minimal interval containing q mass by greedy expansion around mode."""
    prob = np.asarray(prob, float)
    prob = np.clip(prob, 1e-20, None)
    prob = prob / prob.sum()
    idx = int(np.argmax(prob))
    L = R = idx; acc = prob[idx]
    n = prob.size
    while acc < q and (L>0 or R+1<n):
        left = prob[L-1] if L>0 else -1.0
        right= prob[R+1] if R+1<n else -1.0
        if right >= left and R+1<n: R+=1; acc += prob[R]
        elif L>0: L-=1; acc += prob[L]
        else: break
    poc_price = prices[idx]
    va_low = prices[L]; va_high = prices[R]
    return poc_price, va_low, va_high, idx, L, R

# -------------------- MH sampler --------------------
rng = np.random.default_rng(42)

def dirichlet_rw(w, kappa=200.0):
    """Random walk on simplex via Dirichlet centered at w with concentration kappa."""
    alpha = np.maximum(1e-6, w * kappa)
    return rng.dirichlet(alpha)

def mh_fit_mixture(counts, prices, comps, sigma0_bins=2.0, sigma_bounds=(1.0,4.0),
                   iters=1200, burn=400, kappa=200.0, priors=None):
    """
    counts: per-price counts (volume)
    comps: list of per-scale components (each sums ~1)
    priors: dict with 'alpha' (Dirichlet over weights), 'sigma_center' and 'sigma_sd'
    """
    K = len(prices)
    S = len(comps)
    # Normalize comps to probabilities
    comps = [np.clip(c, 1e-20, None) / np.maximum(1e-20, np.sum(c)) for c in comps]
    counts = np.asarray(counts, float)
    N = counts.sum()
    # start
    w = np.ones(S)/S
    sigma = float(np.clip(sigma0_bins, *sigma_bounds))

    # prior setup
    alpha = np.ones(S)
    sigma_center = sigma
    sigma_sd = 0.5
    if priors:
        alpha = np.asarray(priors.get("alpha", alpha), float)
        sigma_center = float(priors.get("sigma_center", sigma_center))
        sigma_sd = float(priors.get("sigma_sd", sigma_sd))

    def build_p(w, sigma):
        mix = np.zeros(K)
        for s in range(S): mix += w[s] * comps[s]
        mix = gaussian_smooth(mix, sigma)
        mix = np.clip(mix, 1e-20, None); mix /= mix.sum()
        return mix

    def logpost(w, sigma):
        p = build_p(w, sigma)
        # multinomial log-likelihood up to constant: sum c_i log p_i
        ll = float(np.sum(counts * np.log(p)))
        # priors: Dirichlet(w|alpha) ~ sum (alpha-1) log w
        lp_w = float(np.sum((alpha - 1.0) * np.log(np.clip(w,1e-20,None))))
        # sigma prior ~ log-normal around sigma_center with sd (rough)
        lp_s = -0.5 * ((np.log(sigma) - np.log(sigma_center))/sigma_sd)**2
        return ll + lp_w + lp_s

    lp = logpost(w, sigma)
    samples_w = []
    samples_sigma = []
    accept_w = accept_s = 0

    for it in range(iters):
        # block 1: weights
        w_prop = dirichlet_rw(w, kappa=kappa)
        lp_prop = logpost(w_prop, sigma)
        # proposal asymmetry correction (Hastings ratio):
        # q(w->w') / q(w'->w) ≈ Dir(w'; kappa*w) / Dir(w; kappa*w')
        def log_dir_pdf(x, alpha):
            return np.sum((alpha-1)*np.log(np.clip(x,1e-20,None))) - (np.sum(np.log(np.arange(1,int(np.sum(alpha)),1))) if False else 0.0)
        # Use approximate symmetric assumption (large kappa) to keep it light:
        if np.log(rng.random()) < (lp_prop - lp):
            w = w_prop; lp = lp_prop; accept_w += 1

        # block 2: sigma (log-rw)
        sig_prop = float(np.clip(np.exp(np.log(sigma) + rng.normal(0, 0.15)), *sigma_bounds))
        lp_prop2 = logpost(w, sig_prop)
        if np.log(rng.random()) < (lp_prop2 - lp):
            sigma = sig_prop; lp = lp_prop2; accept_s += 1

        if it >= burn:
            samples_w.append(w.copy())
            samples_sigma.append(sigma)

    acc_w = accept_w / max(1, iters)
    acc_s = accept_s / max(1, iters)
    W = np.array(samples_w)
    Sg = np.array(samples_sigma)
    # posterior summaries
    w_mean = W.mean(axis=0)
    sig_med = float(np.median(Sg))

    # build posterior draw profiles for POC/VA bands
    draw_ct = min(200, W.shape[0])
    idxs = np.linspace(0, W.shape[0]-1, draw_ct, dtype=int)
    poc_list, val_list, vah_list = [], [], []
    for j in idxs:
        wj = W[j]; sj = Sg[j]
        p = build_p(wj, sj)
        poc, vaL, vaH, _, _, _ = poc_and_va(p, prices, q=0.70)
        poc_list.append(poc); val_list.append(vaL); vah_list.append(vaH)

    post = {
        "w_mean": w_mean.tolist(),
        "sigma_median": sig_med,
        "poc_median": float(np.median(poc_list)),
        "poc_lo50": float(np.percentile(poc_list, 25)),
        "poc_hi50": float(np.percentile(poc_list, 75)),
        "va_low_median": float(np.median(val_list)),
        "va_low50": float(np.percentile(val_list, 25)),
        "va_high50": float(np.percentile(vah_list, 75)),
        "va_high_median": float(np.median(vah_list)),
        "acc_w": acc_w, "acc_s": acc_s
    }
    return post

# -------------------- covariance --------------------
class OnlineCov:
    def __init__(self, hl_seconds=300.0):
        self.lambda_ = math.exp(-math.log(2)/hl_seconds)
        self.mu = None
        self.S = None
    def update(self, x):
        lam = self.lambda_
        x = np.asarray(x, float)
        if self.mu is None:
            self.mu = x.copy()
            self.S = np.eye(x.size)*1e-6
        else:
            self.mu = lam*self.mu + (1-lam)*x
            dx = x - self.mu
            self.S  = lam*self.S + (1-lam)*np.outer(dx, dx)
        # shrinkage
        if _HAVE_SK:
            lw = LedoitWolf().fit(self.S)
            Ssh = lw.covariance_
        else:
            d = np.diag(self.S)
            tau = d.mean()
            alpha = 0.1
            Ssh = (1-alpha)*self.S + alpha*tau*np.eye(self.S.shape[0])
        # λmax
        try:
            lam_max = float(np.linalg.eigvalsh(Ssh).max())
        except Exception:
            lam_max = float(np.nan)
        return lam_max

online_cov = OnlineCov(hl_seconds=300.0)

# -------------------- main loop --------------------
def process_once():
    bins, feat = read_latest_bins_features()
    if bins is None: return None
    prices_ticks = bins["price_ticks"].values.astype(np.int64)
    volumes = bins["volume"].values.astype(float)
    tick_size = float(bins["tick_size"].iloc[0])
    grid_ticks = int(bins["grid_ticks"].iloc[0])
    # Deduplicate by price just in case
    df = pd.DataFrame({"pt":prices_ticks,"v":volumes}).groupby("pt",as_index=False).sum()
    df = df.sort_values("pt")
    counts = df["v"].values
    prices = df["pt"].values * tick_size

    # Build per-price curve (normalized)
    mass = counts.sum()
    if mass <= 0: return None
    curve = counts / mass

    # Components via SWT (or fallback)
    S = 3  # start with 3 scales
    comps = swt_components(curve, wavelet="db4", levels=S)

    # Priors (simple, can be adapted from opening/day-type once you log them)
    priors = {"alpha": np.ones(S), "sigma_center": 2.0, "sigma_sd": 0.5}

    # MH
    post = mh_fit_mixture(
        counts=counts, prices=prices, comps=comps,
        sigma0_bins=float(feat.get("sigma_bins", 2.0)),
        sigma_bounds=(1.0, 4.0), iters=1200, burn=400, kappa=220.0, priors=priors
    )

    # Covariance feature vector (toy: weights + VA width + tails)
    va_width = post["va_high_median"] - post["va_low_median"]
    # crude tails from fallback VA (can be refined using posterior draws)
    tail_up = max(0.0, prices.max() - post["va_high_median"])
    tail_dn = max(0.0, post["va_low_median"] - prices.min())
    x = np.array([*post["w_mean"], va_width, tail_up, tail_dn], float)
    lam_max = online_cov.update(x)

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
        # pad weights to 4 cols for Sierra reader (w4=0 if S=3)
        "weights_extended": (post["w_mean"] + [0.0])[:4]
    }
    write_results(res)
    return res

def main():
    print("Bayes Sidecar: watching", DATA_DIR)
    last_sig = None
    while True:
        try:
            sig = None
            for f in ("feed_bins.csv","feed_features.csv"):
                p = os.path.join(DATA_DIR, f)
                if os.path.exists(p):
                    t = os.path.getmtime(p)
                    sig = t if sig is None else max(sig, t)
            if sig and sig != last_sig:
                out = process_once()
                last_sig = sig
        except KeyboardInterrupt:
            break
        except Exception as e:
            print("Error:", e)
        time.sleep(1.0)

if __name__ == "__main__":
    main()
