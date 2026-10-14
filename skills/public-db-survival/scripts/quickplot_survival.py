#!/usr/bin/env python3
"""Dependency-light (matplotlib + numpy only) survival quicklook generator.

Produces the three Required Figures from references/cns-method-guide.md:
  - expression box/violin plot with sample counts
  - KM curve with risk table, censoring ticks, group sizes, HR/CI and
    log-rank p value
  - Cox proportional-hazards forest plot (univariate per covariate, plus a
    multivariate model when more than one covariate is supplied)

The Cox model is fit with a from-scratch Newton-Raphson solver on the
Breslow-tie-approximation partial likelihood (numpy only, no scipy/lifelines).
This is an approximation intended for quicklook/exploratory use; for a
formal Efron-tie-corrected model use scripts/run_survival.R /
scripts/run_cox_forest.R when R + survival (+ survminer/ggplot2) are
available.
"""
import argparse
import csv
import json
import math
import statistics
from pathlib import Path

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as exc:
    raise SystemExit(
        "quickplot_survival.py requires matplotlib to export PNG/PDF. "
        "Install it with: python3 -m pip install matplotlib"
    ) from exc

try:
    import numpy as np
except Exception as exc:
    raise SystemExit(
        "quickplot_survival.py requires numpy for the Cox model fit. "
        "Install it with: python3 -m pip install numpy"
    ) from exc


TEXT = "#173B5F"
GRID = "#DCEAF5"
HIGH = "#E45B5B"
LOW = "#4EA3F1"
NEUTRAL = "#667085"


def apply_style():
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 16,
            "axes.titlesize": 24,
            "axes.labelsize": 21,
            "xtick.labelsize": 16,
            "ytick.labelsize": 16,
            "legend.fontsize": 15,
            "figure.titlesize": 22,
            "axes.edgecolor": "#B9D7EE",
            "axes.labelcolor": TEXT,
            "xtick.color": TEXT,
            "ytick.color": TEXT,
            "text.color": TEXT,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save_plot(fig, path_base):
    path_base = Path(path_base)
    path_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path_base.with_suffix(".png"), dpi=320, bbox_inches="tight", facecolor="white")
    fig.savefig(path_base.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def read_csv(path):
    return list(csv.DictReader(Path(path).read_text(encoding="utf-8-sig").splitlines()))


def coerce_event(value):
    """Map 0/1, alive/dead, censored/event style values to an int 0/1."""
    s = str(value).strip().lower()
    if s in ("1", "dead", "deceased", "event", "yes", "true"):
        return 1
    if s in ("0", "alive", "censored", "no", "false"):
        return 0
    return int(float(value))


# ---------------------------------------------------------------------------
# Kaplan-Meier + log-rank
# ---------------------------------------------------------------------------

def km_curve(rows, time_col, status_col):
    """Return (times, surv, censor_times, censor_surv, risk_table_fn).

    times/surv describe the right-continuous KM step function starting at
    (0, 1.0). censor_times/censor_surv give the tick-mark positions for
    censored observations. risk_table_fn(t) returns the number still at
    risk (not yet failed/censored) at time t.
    """
    data = sorted((float(r[time_col]), coerce_event(r[status_col])) for r in rows)
    n = len(data)
    surv = 1.0
    pts = [(0.0, 1.0)]
    censor_times, censor_surv = [], []
    at_risk = n
    for t, event in data:
        if event == 1 and at_risk > 0:
            surv *= (at_risk - 1) / at_risk
            pts.append((t, surv))
        else:
            censor_times.append(t)
            censor_surv.append(surv)
        at_risk -= 1

    def risk_at(t):
        return sum(1 for (tt, _) in data if tt >= t)

    return pts, censor_times, censor_surv, risk_at, data


def group_survival_summary(rows, group_col, time_col, status_col):
    """Per-group n, events, median survival ('not reached' if the KM curve
    never crosses 0.5, per references/pitfalls.md #15), and last follow-up."""
    summary = {}
    for group in sorted({r[group_col] for r in rows}):
        group_rows = [r for r in rows if r[group_col] == group]
        pts, _, _, _, data = km_curve(group_rows, time_col, status_col)
        last_time = data[-1][0] if data else float("nan")
        median = next((t for t, s in pts if s <= 0.5), None)
        summary[group] = {
            "n": len(group_rows),
            "events": sum(coerce_event(r[status_col]) for r in group_rows),
            "median_survival": "not reached" if median is None else round(median, 4),
            "last_followup": round(last_time, 4),
        }
    return summary


def logrank_test(rows, group_col, time_col, status_col, groups):
    """Two-group Mantel-Haenszel log-rank test (no scipy).

    Returns (chisq, pvalue, observed) for the first listed group ("O1").
    """
    if len(groups) != 2:
        return float("nan"), float("nan"), None
    g1, g2 = groups
    records = []
    for r in rows:
        if r[group_col] not in (g1, g2):
            continue
        records.append((float(r[time_col]), coerce_event(r[status_col]), r[group_col] == g1))
    event_times = sorted({t for (t, e, _) in records if e == 1})
    o1 = e1 = var = 0.0
    for t in event_times:
        at_risk = [rec for rec in records if rec[0] >= t]
        n_j = len(at_risk)
        n1_j = sum(1 for rec in at_risk if rec[2])
        deaths = [rec for rec in at_risk if rec[0] == t and rec[1] == 1]
        d_j = len(deaths)
        d1_j = sum(1 for rec in deaths if rec[2])
        if n_j <= 1 or d_j == 0:
            continue
        o1 += d1_j
        e1 += d_j * n1_j / n_j
        if n_j > 1:
            var += d_j * (n1_j / n_j) * (1 - n1_j / n_j) * (n_j - d_j) / (n_j - 1)
    if var <= 0:
        return float("nan"), float("nan"), o1
    chisq = (o1 - e1) ** 2 / var
    pvalue = 1.0 - math.erf(math.sqrt(chisq / 2.0))
    return chisq, pvalue, o1


# ---------------------------------------------------------------------------
# Cox proportional hazards (Newton-Raphson, Breslow tie approximation)
# ---------------------------------------------------------------------------

def _norm_cdf(x):
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def coxph_fit(times, events, X, max_iter=50, tol=1e-6):
    """Fit a Cox PH model via Newton-Raphson on the Breslow partial likelihood.

    times: (n,) array of follow-up times.
    events: (n,) array of 0/1 event indicators.
    X: (n, p) covariate matrix (already numeric, no intercept column).

    Returns a dict with coef, se, hr, lower95, upper95, pvalue (each length p),
    n, n_events, converged. Ties are handled with the Breslow approximation
    (simpler than R's default Efron correction; adequate for a quicklook
    figure, not a substitute for scripts/run_cox_forest.R).
    """
    times = np.asarray(times, dtype=float)
    events = np.asarray(events, dtype=float)
    X = np.asarray(X, dtype=float)
    n, p = X.shape

    order = np.argsort(times)
    t_sorted = times[order]
    e_sorted = events[order]
    X_sorted = X[order]

    beta = np.zeros(p)
    converged = False
    for _ in range(max_iter):
        eta = X_sorted @ beta
        eta -= eta.max()  # numerical stability, cancels in ratios
        w = np.exp(eta)
        # Suffix sums: risk set for row i is all rows with time >= t_sorted[i],
        # which (after ascending sort) is exactly the suffix i..n-1.
        S0 = np.cumsum(w[::-1])[::-1]
        S1 = np.cumsum((w[:, None] * X_sorted)[::-1], axis=0)[::-1]
        outer = np.einsum("ni,nj->nij", X_sorted, X_sorted) * w[:, None, None]
        S2 = np.cumsum(outer[::-1], axis=0)[::-1]

        grad = np.zeros(p)
        hess = np.zeros((p, p))
        event_idx = np.where(e_sorted == 1)[0]
        for i in event_idx:
            s1_over_s0 = S1[i] / S0[i]
            grad += X_sorted[i] - s1_over_s0
            hess -= S2[i] / S0[i] - np.outer(s1_over_s0, s1_over_s0)

        try:
            step = np.linalg.solve(hess, grad)
        except np.linalg.LinAlgError:
            step = np.linalg.lstsq(hess, grad, rcond=None)[0]
        beta_new = beta - step
        if np.max(np.abs(beta_new - beta)) < tol:
            beta = beta_new
            converged = True
            break
        beta = beta_new

    # Recompute Hessian at the final beta for inference.
    eta = X_sorted @ beta
    eta -= eta.max()
    w = np.exp(eta)
    S0 = np.cumsum(w[::-1])[::-1]
    S1 = np.cumsum((w[:, None] * X_sorted)[::-1], axis=0)[::-1]
    outer = np.einsum("ni,nj->nij", X_sorted, X_sorted) * w[:, None, None]
    S2 = np.cumsum(outer[::-1], axis=0)[::-1]
    hess = np.zeros((p, p))
    for i in np.where(e_sorted == 1)[0]:
        s1_over_s0 = S1[i] / S0[i]
        hess -= S2[i] / S0[i] - np.outer(s1_over_s0, s1_over_s0)

    fisher_info = -hess
    try:
        cov = np.linalg.inv(fisher_info)
    except np.linalg.LinAlgError:
        cov = np.linalg.pinv(fisher_info)
    se = np.sqrt(np.clip(np.diag(cov), 0, None))
    z = np.divide(beta, se, out=np.full_like(beta, np.nan), where=se > 0)
    pvalue = np.array([2 * (1 - _norm_cdf(abs(zi))) if not math.isnan(zi) else float("nan") for zi in z])
    hr = np.exp(beta)
    lower95 = np.exp(beta - 1.9599639845400545 * se)
    upper95 = np.exp(beta + 1.9599639845400545 * se)
    return {
        "coef": beta,
        "se": se,
        "hr": hr,
        "lower95": lower95,
        "upper95": upper95,
        "pvalue": pvalue,
        "n": n,
        "n_events": int(e_sorted.sum()),
        "converged": converged,
    }


def schoenfeld_ph_test(times, events, X, beta, var_names):
    """Approximate proportional-hazards trend test (numpy only, no scipy).

    Computes Schoenfeld residuals r_i = x_i - S1(t_i)/S0(t_i) at each event
    time from the already-fitted Cox model, then tests the Pearson
    correlation between each covariate's residual series and the rank of
    event time via a Fisher z-transform (normal approximation). This is the
    same "does the covariate effect drift over time" question that
    survival::cox.zph (Grambsch & Therneau 1994) answers, but WITHOUT the
    score-covariance variance-scaling used there, so treat this as an
    exploratory proxy, not a formal PH test; prefer scripts/run_survival.R /
    scripts/run_cox_forest.R (survival::cox.zph) for a formal result.
    """
    times = np.asarray(times, dtype=float)
    events = np.asarray(events, dtype=float)
    X = np.asarray(X, dtype=float)
    order = np.argsort(times)
    t_sorted = times[order]
    e_sorted = events[order]
    X_sorted = X[order]

    eta = X_sorted @ beta
    eta -= eta.max()
    w = np.exp(eta)
    S0 = np.cumsum(w[::-1])[::-1]
    S1 = np.cumsum((w[:, None] * X_sorted)[::-1], axis=0)[::-1]

    event_idx = np.where(e_sorted == 1)[0]
    n_events = len(event_idx)
    results = []
    if n_events < 4:
        for name in var_names:
            results.append({"variable": name, "correlation": float("nan"), "pvalue": float("nan"), "n_events": n_events})
        return results

    resid = X_sorted[event_idx] - S1[event_idx] / S0[event_idx][:, None]
    etimes = t_sorted[event_idx]
    ranks = np.argsort(np.argsort(etimes)).astype(float) + 1.0

    for j, name in enumerate(var_names):
        r = resid[:, j]
        if np.std(r) == 0 or np.std(ranks) == 0:
            results.append({"variable": name, "correlation": 0.0, "pvalue": float("nan"), "n_events": n_events})
            continue
        corr = float(np.corrcoef(ranks, r)[0, 1])
        corr = max(min(corr, 0.999999), -0.999999)
        z = math.atanh(corr) * math.sqrt(max(n_events - 3, 1))
        pvalue = 2 * (1 - _norm_cdf(abs(z)))
        results.append({"variable": name, "correlation": corr, "pvalue": pvalue, "n_events": n_events})
    return results


def ph_status_label(ph_rows):
    """Turn schoenfeld_ph_test() rows into a plain-language PH status string."""
    valid = [r for r in ph_rows if not math.isnan(r["pvalue"])]
    if not valid:
        return "not assessed (too few events for a trend test)", False, []
    violated = [r["variable"] for r in valid if r["pvalue"] < 0.05]
    if violated:
        return (
            f"VIOLATED for covariate(s) {', '.join(violated)} (approximate rank-correlation test p<0.05): "
            "HR(s) should be interpreted with caution as time-averaged approximations.",
            True,
            violated,
        )
    return (
        "satisfied for " + ", ".join(r["variable"] for r in valid) + " (approximate rank-correlation test p>=0.05)",
        False,
        [],
    )


def build_cox_table(rows, time_col, status_col, covariates):
    """Fit univariate Cox per covariate, plus a joint multivariate model
    when more than one covariate is given. Returns a list of dict rows
    (model, variable, HR, lower95, upper95, pvalue, n, n_events)."""
    times = np.array([float(r[time_col]) for r in rows])
    events = np.array([coerce_event(r[status_col]) for r in rows])
    if events.sum() == 0:
        raise ValueError(
            "zero events observed in this cohort (all-censored data); a Cox model has no information to fit "
            "and any reported hazard ratio would be an artifact, not a result"
        )
    table = []
    for cov in covariates:
        x = np.array([[float(r[cov])] for r in rows])
        fit = coxph_fit(times, events, x)
        table.append(
            {
                "model": "univariate",
                "variable": cov,
                "HR": fit["hr"][0],
                "lower95": fit["lower95"][0],
                "upper95": fit["upper95"][0],
                "pvalue": fit["pvalue"][0],
                "n": fit["n"],
                "n_events": fit["n_events"],
                "converged": fit["converged"],
            }
        )
    if len(covariates) > 1:
        X = np.array([[float(r[c]) for c in covariates] for r in rows])
        fit = coxph_fit(times, events, X)
        for j, cov in enumerate(covariates):
            table.append(
                {
                    "model": "multivariate",
                    "variable": cov,
                    "HR": fit["hr"][j],
                    "lower95": fit["lower95"][j],
                    "upper95": fit["upper95"][j],
                    "pvalue": fit["pvalue"][j],
                    "n": fit["n"],
                    "n_events": fit["n_events"],
                    "converged": fit["converged"],
                }
            )
    return table


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def write_km(rows, group_col, time_col, status_col, out_base, cutoff_note=""):
    apply_style()
    groups = sorted({r[group_col] for r in rows})
    colors = {"Low": LOW, "High": HIGH}
    group_n = {g: sum(1 for r in rows if r[group_col] == g) for g in groups}

    hr_label = ""
    logrank_label = ""
    hr_val = hr_lower = hr_upper = None
    pvalue = float("nan")
    total_events = sum(coerce_event(r[status_col]) for r in rows)
    if len(groups) == 2 and total_events == 0:
        hr_label = "not estimable (zero events observed in this cohort; all-censored data carries no survival information)"
    elif len(groups) == 2:
        chisq, pvalue, _ = logrank_test(rows, group_col, time_col, status_col, ["Low", "High"] if set(groups) == {"Low", "High"} else groups)
        if not math.isnan(pvalue):
            logrank_label = f"Log-rank p = {pvalue:.3g}"
        try:
            x = np.array([[1.0 if r[group_col] == "High" else 0.0] for r in rows])
            times = np.array([float(r[time_col]) for r in rows])
            events = np.array([coerce_event(r[status_col]) for r in rows])
            fit = coxph_fit(times, events, x)
            hr_val, hr_lower, hr_upper = float(fit["hr"][0]), float(fit["lower95"][0]), float(fit["upper95"][0])
            hr_label = f"HR(High vs Low) = {hr_val:.2f} (95% CI {hr_lower:.2f}-{hr_upper:.2f})"
        except Exception:
            hr_label = ""

    time_max = max(float(r[time_col]) for r in rows)
    breaks = np.linspace(0, time_max, 6)

    fig, (ax, ax_table) = plt.subplots(
        2, 1, figsize=(8.8, 7.8), gridspec_kw={"height_ratios": [3, 1]}, sharex=False
    )
    ax.set_facecolor("#FBFDFF")
    risk_fns = {}
    for group in groups:
        group_rows = [r for r in rows if r[group_col] == group]
        pts, cens_t, cens_s, risk_at, _ = km_curve(group_rows, time_col, status_col)
        x_vals = [p[0] for p in pts]
        y_vals = [p[1] for p in pts]
        color = colors.get(group, NEUTRAL)
        ax.step(x_vals, y_vals, where="post", linewidth=3.0, color=color, label=f"{group} n={group_n[group]}")
        if cens_t:
            ax.plot(cens_t, cens_s, marker="+", linestyle="none", markersize=9, markeredgewidth=1.6, color=color)
        risk_fns[group] = risk_at

    ax.set_title("Kaplan-Meier curve", pad=18, loc="center", fontweight="bold")
    subtitle = " | ".join(x for x in [cutoff_note, hr_label, logrank_label] if x)
    if subtitle:
        ax.text(0.5, 1.02, subtitle, transform=ax.transAxes, ha="center", va="bottom", fontsize=11.5, color=TEXT)
    ax.set_xlabel("")
    ax.set_ylabel("Survival probability", labelpad=12)
    ax.set_ylim(0, 1.05)
    ax.set_xlim(0, time_max * 1.02)
    ax.grid(True, color=GRID, linewidth=0.9)
    ax.legend(frameon=False, loc="upper right")

    # Extend the left edge of the risk-table x-axis so the at-risk count
    # centered at t=0 has room to render without colliding with the group
    # name tick label to its left (both otherwise sit right at x=0).
    table_xlim = (ax.get_xlim()[0] - time_max * 0.09, ax.get_xlim()[1])
    ax_table.set_xlim(table_xlim)
    ax_table.set_ylim(-0.5, len(groups) - 0.5)
    ax_table.set_yticks(range(len(groups)))
    ax_table.set_yticklabels(list(reversed(groups)))
    ax_table.set_xlabel("Time", labelpad=10)
    for spine in ("top", "right", "left"):
        ax_table.spines[spine].set_visible(False)
    for i, group in enumerate(reversed(groups)):
        for t in breaks:
            n_at_risk = risk_fns[group](t)
            ax_table.text(t, i, str(n_at_risk), ha="center", va="center", fontsize=12, color=colors.get(group, NEUTRAL))
    ax_table.set_title("Number at risk", loc="left", fontsize=13, pad=6)
    ax_table.tick_params(left=False, pad=12)

    fig.tight_layout()
    save_plot(fig, out_base)
    return {
        "logrank_p": pvalue if len(groups) == 2 else None,
        "hr_label": hr_label,
        "hr": hr_val,
        "hr_lower95": hr_lower,
        "hr_upper95": hr_upper,
    }


def write_box(rows, expr_col, type_col, out_base):
    apply_style()
    groups = sorted({r[type_col] for r in rows})
    data = [[float(r[expr_col]) for r in rows if r[type_col] == group] for group in groups]
    colors = [LOW, HIGH, "#4CB984", "#8B72D6"]
    fig, ax = plt.subplots(figsize=(7.8, 6.4))
    box = ax.boxplot(data, patch_artist=True, widths=0.55, medianprops={"color": TEXT, "linewidth": 2.2})
    ax.set_xticks(range(1, len(groups) + 1))
    ax.set_xticklabels([f"{g}\nn={len(vals)}" for g, vals in zip(groups, data)])
    for patch, color in zip(box["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.24)
        patch.set_edgecolor(color)
        patch.set_linewidth(1.8)
    for whisker in box["whiskers"]:
        whisker.set_color("#7BAED5")
        whisker.set_linewidth(1.5)
    for cap in box["caps"]:
        cap.set_color("#7BAED5")
        cap.set_linewidth(1.5)
    for i, vals in enumerate(data, start=1):
        offsets = [((idx % 7) - 3) * 0.018 for idx in range(len(vals))]
        ax.scatter([i + off for off in offsets], vals, s=22, color=colors[(i - 1) % len(colors)], alpha=0.55, linewidths=0)
    ax.set_title(f"{expr_col} expression by {type_col}", pad=18, loc="center", fontweight="bold")
    ax.set_xlabel(type_col, labelpad=12)
    ax.set_ylabel("Expression", labelpad=12)
    ax.grid(True, axis="y", color=GRID, linewidth=0.9)
    fig.tight_layout()
    save_plot(fig, out_base)


def write_cox_forest(table, out_base):
    """Draw a forest plot (matplotlib) from build_cox_table() rows."""
    apply_style()
    variables = list(dict.fromkeys(row["variable"] for row in table))
    n = len(variables)
    has_multi = any(row["model"] == "multivariate" for row in table)
    model_color = {"univariate": NEUTRAL, "multivariate": HIGH}
    model_offset = {"univariate": -0.15, "multivariate": 0.15} if has_multi else {"univariate": 0.0}

    fig, ax = plt.subplots(figsize=(11.4, max(3.4, 1.0 + 0.7 * n)))
    y_base = {var: n - i for i, var in enumerate(variables)}

    # Reserve the right ~45% of the log-scale x-axis for the HR/CI/p text so
    # labels never collide with the point estimates or the legend.
    all_bounds = [row["lower95"] for row in table] + [row["upper95"] for row in table] + [1.0]
    log_bounds = [math.log10(v) for v in all_bounds if v > 0]
    lo, hi = min(log_bounds), max(log_bounds)
    span = max(hi - lo, 0.05)
    xlim = (10 ** (lo - 0.15 * span), 10 ** (hi + 1.15 * span))

    for row in table:
        y = y_base[row["variable"]] + model_offset[row["model"]]
        color = model_color[row["model"]]
        ax.plot([row["lower95"], row["upper95"]], [y, y], color=color, linewidth=2.2, solid_capstyle="round")
        ax.plot(row["HR"], y, marker="D", markersize=7, color=color)
        label = f"{row['HR']:.2f} ({row['lower95']:.2f}-{row['upper95']:.2f}), p={row['pvalue']:.3g}"
        ax.annotate(label, xy=(10 ** (hi + 0.2 * span), y), xytext=(0, 0), textcoords="offset points",
                    va="center", ha="left", fontsize=11, color=color)

    ax.axvline(1, linestyle="--", color="#A7C6DD", linewidth=1.2)
    ax.set_xscale("log")
    ax.set_xlim(xlim)
    ax.set_yticks(list(y_base.values()))
    ax.set_yticklabels(list(y_base.keys()))
    ax.set_ylim(0.3, n + 0.7)
    ax.set_xlabel("Hazard ratio (log scale)", labelpad=12)
    ax.set_title("Cox model forest plot", pad=18, loc="center", fontweight="bold")
    if has_multi:
        handles = [
            plt.Line2D([0], [0], color=NEUTRAL, marker="D", linestyle="", label="Univariate"),
            plt.Line2D([0], [0], color=HIGH, marker="D", linestyle="", label="Multivariate"),
        ]
        ax.legend(handles=handles, frameon=False, loc="lower right")
    ax.grid(True, axis="x", color=GRID, linewidth=0.9)
    fig.tight_layout()
    save_plot(fig, out_base)


def main():
    parser = argparse.ArgumentParser(description="Survival quick plots with PNG/PDF export.")
    parser.add_argument("--table", required=True)
    parser.add_argument("--gene", default="GeneA")
    parser.add_argument("--time-col", default="survival_month")
    parser.add_argument("--status-col", default="status")
    parser.add_argument("--type-col", default="type")
    parser.add_argument("--cutoff-method", default="median", choices=["median", "tertile", "quartile"])
    parser.add_argument("--covariates", default="", help="Comma-separated extra numeric covariate columns for the Cox forest plot (target gene is always included)")
    parser.add_argument("--out", default="output")
    args = parser.parse_args()
    out = Path(args.out)
    (out / "figures").mkdir(parents=True, exist_ok=True)
    (out / "results").mkdir(parents=True, exist_ok=True)
    (out / "run").mkdir(parents=True, exist_ok=True)
    rows = read_csv(args.table)
    if not rows:
        raise SystemExit(f"Input table '{args.table}' has no data rows.")

    required_cols = [args.gene, args.time_col, args.status_col, args.type_col]
    missing_cols = [c for c in required_cols if c not in rows[0]]
    if missing_cols:
        raise SystemExit(
            f"Column(s) not found in '{args.table}': {', '.join(missing_cols)}. "
            f"Available columns: {', '.join(rows[0].keys())}"
        )

    n_before_numeric = len(rows)
    bad_rows = []
    clean_rows = []
    for r in rows:
        try:
            float(r[args.gene])
            float(r[args.time_col])
            coerce_event(r[args.status_col])
            clean_rows.append(r)
        except (ValueError, TypeError):
            bad_rows.append(r)
    if bad_rows:
        print(
            f"Note: dropped {len(bad_rows)} row(s) with missing/non-numeric '{args.gene}', "
            f"'{args.time_col}', or unmappable '{args.status_col}' values "
            f"(kept {len(clean_rows)}/{n_before_numeric})."
        )
    rows = clean_rows
    if not rows:
        raise SystemExit(
            f"No usable rows remain after removing missing/non-numeric '{args.gene}'/'{args.time_col}'/"
            f"'{args.status_col}' values."
        )

    values = [float(r[args.gene]) for r in rows]
    method = args.cutoff_method
    if method == "median":
        cutoff = statistics.median(values)
        for r in rows:
            r["expression_group"] = "High" if float(r[args.gene]) >= cutoff else "Low"
        cutoff_note = f"Median split at {cutoff:.4g}"
    else:
        probs = (1 / 3, 2 / 3) if method == "tertile" else (0.25, 0.75)
        sv = sorted(values)

        def quantile(p):
            idx = p * (len(sv) - 1)
            lo, hi = int(math.floor(idx)), int(math.ceil(idx))
            if lo == hi:
                return sv[lo]
            frac = idx - lo
            return sv[lo] * (1 - frac) + sv[hi] * frac

        q_lo, q_hi = quantile(probs[0]), quantile(probs[1])
        for r in rows:
            v = float(r[args.gene])
            r["expression_group"] = "Low" if v <= q_lo else ("High" if v >= q_hi else "Middle")
        n_before = len(rows)
        rows = [r for r in rows if r["expression_group"] != "Middle"]
        cutoff = (q_lo, q_hi)
        cutoff_note = f"{method.title()} split (dropped {n_before - len(rows)} middle samples)"

    km_stats = write_km(rows, "expression_group", args.time_col, args.status_col, out / "figures" / "km_curve", cutoff_note)
    write_box(rows, args.gene, args.type_col, out / "figures" / "expression_boxplot")
    group_km_summary = group_survival_summary(rows, "expression_group", args.time_col, args.status_col)

    covariates = [args.gene] + [c.strip() for c in args.covariates.split(",") if c.strip()]
    cox_table = []
    try:
        cox_table = build_cox_table(rows, args.time_col, args.status_col, covariates)
        write_cox_forest(cox_table, out / "figures" / "cox_forest")
    except Exception as exc:
        (out / "run" / "cox_forest_error.txt").write_text(f"Cox forest plot failed: {exc}\n", encoding="utf-8")
        print(f"Warning: Cox forest plot failed ({exc}); see output/run/cox_forest_error.txt")

    # Proportional-hazards check (references/pitfalls.md #2): the R route uses
    # survival::cox.zph; this Python route uses the from-scratch approximate
    # Schoenfeld-residual rank-correlation test defined above so a hazard
    # ratio is never reported without a PH check, even without R available.
    ph_rows = []
    ph_group_label = ph_cont_label = ph_multi_label = None
    ph_group_violated = ph_cont_violated = ph_multi_violated = False
    epv = None
    try:
        times_arr = np.array([float(r[args.time_col]) for r in rows])
        events_arr = np.array([coerce_event(r[args.status_col]) for r in rows])

        x_group = np.array([[1.0 if r["expression_group"] == "High" else 0.0] for r in rows])
        group_fit = coxph_fit(times_arr, events_arr, x_group)
        group_ph = schoenfeld_ph_test(times_arr, events_arr, x_group, group_fit["coef"], ["expression_group(High vs Low)"])
        ph_group_label, ph_group_violated, _ = ph_status_label(group_ph)
        ph_rows += [dict(model="expression_group (categorical)", **r) for r in group_ph]

        x_cont = np.array([[float(r[args.gene])] for r in rows])
        cont_fit = coxph_fit(times_arr, events_arr, x_cont)
        cont_ph = schoenfeld_ph_test(times_arr, events_arr, x_cont, cont_fit["coef"], [args.gene])
        ph_cont_label, ph_cont_violated, _ = ph_status_label(cont_ph)
        ph_rows += [dict(model=f"{args.gene} (continuous)", **r) for r in cont_ph]

        if len(covariates) > 1:
            X_multi = np.array([[float(r[c]) for c in covariates] for r in rows])
            multi_fit = coxph_fit(times_arr, events_arr, X_multi)
            multi_ph = schoenfeld_ph_test(times_arr, events_arr, X_multi, multi_fit["coef"], covariates)
            ph_multi_label, ph_multi_violated, _ = ph_status_label(multi_ph)
            ph_rows += [dict(model="multivariate", **r) for r in multi_ph]
            epv = multi_fit["n_events"] / len(covariates)
    except Exception as exc:
        (out / "run" / "ph_test_error.txt").write_text(f"PH check failed: {exc}\n", encoding="utf-8")
        print(f"Warning: PH check failed ({exc}); see output/run/ph_test_error.txt")

    if ph_rows:
        with (out / "results" / "ph_test.tsv").open("w", newline="", encoding="utf-8") as handle:
            fieldnames = ["model", "variable", "correlation", "pvalue", "n_events"]
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for row in ph_rows:
                writer.writerow(row)

    with (out / "results" / "expression_group.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    group_n = {g: sum(1 for r in rows if r["expression_group"] == g) for g in sorted({r["expression_group"] for r in rows})}
    group_events = {
        g: sum(coerce_event(r[args.status_col]) for r in rows if r["expression_group"] == g) for g in group_n
    }
    with (out / "results" / "survival_summary.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["item", "value"])
        writer.writerow(["gene", args.gene])
        writer.writerow(["cutoff_method", method])
        writer.writerow(["cutoff_note", cutoff_note])
        writer.writerow(["samples", len(rows)])
        for g in group_n:
            writer.writerow([f"n_{g}", group_n[g]])
            writer.writerow([f"events_{g}", group_events[g]])
            writer.writerow([f"median_survival_{g}", group_km_summary[g]["median_survival"]])
            writer.writerow([f"last_followup_{g}", group_km_summary[g]["last_followup"]])
        writer.writerow(["logrank_p", km_stats.get("logrank_p")])
        writer.writerow(["hazard_ratio", km_stats.get("hr_label")])
        writer.writerow(["ph_status_group", ph_group_label])
        writer.writerow(["ph_status_continuous", ph_cont_label])
        if epv is not None:
            writer.writerow(["epv_multivariate", round(epv, 2)])
            writer.writerow(["ph_status_multivariate", ph_multi_label])

    if cox_table:
        with (out / "results" / "cox_forest.tsv").open("w", newline="", encoding="utf-8") as handle:
            fieldnames = ["model", "variable", "HR", "lower95", "upper95", "pvalue", "n", "n_events", "converged"]
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for row in cox_table:
                writer.writerow(row)

    cutoff_repr = cutoff if not isinstance(cutoff, tuple) else list(cutoff)
    (out / "summary.tsv").write_text(
        f"item\tvalue\nsamples\t{len(rows)}\ncutoff_method\t{method}\ncutoff\t{cutoff_repr}\n", encoding="utf-8"
    )
    (out / "run" / "qc.json").write_text(
        json.dumps({"samples": len(rows), "cutoff_method": method, "cutoff": cutoff_repr, "quicklook": True}, indent=2),
        encoding="utf-8",
    )

    # --- plain-language conclusion + pitfalls-linked caveats (only what applies) ---
    logrank_p = km_stats.get("logrank_p")
    median_lines = "; ".join(
        f"{g} n={group_km_summary[g]['n']} (events={group_km_summary[g]['events']}), "
        f"median survival={group_km_summary[g]['median_survival']} (last follow-up {group_km_summary[g]['last_followup']})"
        for g in group_km_summary
    )

    caveats = []
    if ph_group_violated or ph_cont_violated:
        caveats.append(
            "the approximate proportional-hazards check flagged a time trend for at least one model above; "
            "read the HR as a time-averaged approximation rather than a constant effect (references/pitfalls.md #2)"
        )
    if epv is not None and epv < 10:
        caveats.append(
            f"events-per-variable for the multivariate model is {epv:.1f}, below the ~10 EPV heuristic "
            "(Peduzzi et al. 1996), so the multivariate coefficients/CIs may be unstable (references/pitfalls.md #8)"
        )
    if any(group_km_summary[g]["median_survival"] == "not reached" for g in group_km_summary):
        caveats.append(
            "median survival was not reached in at least one group within observed follow-up, so treat any "
            "such median as an extrapolation beyond the data (references/pitfalls.md #15)"
        )
    caveats.append(
        "this Cox fit uses a from-scratch Breslow-tie-approximation Newton-Raphson solver and an approximate "
        "(non-variance-scaled) proportional-hazards proxy test, both exploratory; prefer scripts/run_survival.R / "
        "scripts/run_cox_forest.R (survival::coxph + survival::cox.zph) for a formal report"
    )
    caveat_txt = " Caveats: " + "; ".join(caveats) + "." if caveats else ""

    if hr_val := km_stats.get("hr"):
        direction = "worse" if hr_val > 1 else "better"
        sig_txt = (
            f"a statistically significant difference (log-rank p={logrank_p:.3g})"
            if logrank_p is not None and not math.isnan(logrank_p) and logrank_p < 0.05
            else f"no statistically significant difference (log-rank p={logrank_p:.3g})"
            if logrank_p is not None and not math.isnan(logrank_p)
            else "an indeterminate log-rank result"
        )
        conclusion = (
            f"In this cohort, patients with High {args.gene} expression showed {direction} survival than "
            f"patients with Low expression, with {sig_txt} and a Cox hazard ratio of {hr_val:.2f} "
            f"(95% CI {km_stats.get('hr_lower95'):.2f}-{km_stats.get('hr_upper95'):.2f}) for High vs Low."
            f"{caveat_txt} This is a single-cohort, retrospective association and does not by itself establish "
            "causality, clinical utility, or independence from clinical confounders (see references/pitfalls.md)."
        )
    else:
        conclusion = (
            "A Cox hazard ratio could not be estimated for this run (see output/run/cox_forest_error.txt or "
            "output/run/ph_test_error.txt if present); review sample size, event count, and group balance before "
            "drawing any conclusion."
        )

    report_lines = [
        "# Public database survival quicklook",
        "",
        f"Target gene: {args.gene}",
        f"Cutoff rule: {cutoff_note} ({method} is a pre-specified rule; low cutoff-selection-bias risk per references/pitfalls.md #1)",
        f"Samples used: {len(rows)}",
        f"Group sizes and median survival: {median_lines}",
        f"Log-rank test: p={logrank_p:.3g}" if logrank_p is not None and not math.isnan(logrank_p) else "Log-rank test: not computed",
        f"Cox HR (High vs Low): {km_stats.get('hr_label') or 'not computed'}",
        f"Proportional-hazards check (expression_group model, approximate): PH assumption {ph_group_label or 'not assessed'}",
        f"Proportional-hazards check (continuous {args.gene} model, approximate): PH assumption {ph_cont_label or 'not assessed'}",
    ]
    if epv is not None:
        report_lines.append(
            f"Events-per-variable (multivariate model, {len(covariates)} covariates): {epv:.1f}"
            + (" -- BELOW the ~10 EPV heuristic (Peduzzi et al. 1996)" if epv < 10 else " -- meets the ~10 EPV heuristic")
        )
        report_lines.append(
            f"Proportional-hazards check (multivariate model, approximate): PH assumption {ph_multi_label or 'not assessed'}"
        )
    report_lines += [
        "",
        f"Generated {method}-cutoff KM curve (with risk table and censoring ticks), expression boxplot, "
        f"and a Cox forest plot for `{args.gene}`"
        + (f" plus covariates: {', '.join(covariates[1:])}" if len(covariates) > 1 else "")
        + ". Figures are exported as PNG and PDF.",
        "",
        "## Conclusion",
        "",
        conclusion,
        "",
        "## Limitations",
        "",
        "- This quicklook Cox model uses a Breslow-tie approximation fit with numpy only (no scipy/lifelines); "
        "treat it as exploratory and prefer scripts/run_cox_forest.R for a formal Efron-tie-corrected model when R is available.",
        "- The proportional-hazards check above is an approximate Schoenfeld-residual rank-correlation proxy, not the "
        "formal survival::cox.zph (Grambsch & Therneau) test; prefer the R route for a formal PH check.",
        "- Mark small or simulated cohorts as exploratory; this KM/log-rank/Cox result does not by itself establish clinical utility.",
        "- Read references/pitfalls.md before reporting any hazard ratio or prognostic claim (cutoff-selection bias, "
        "unchecked proportional hazards, immortal time bias, competing risks, and more).",
    ]
    (out / "report.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    print(f"wrote survival quicklook outputs to {out}")


if __name__ == "__main__":
    main()
