import pandas as pd
import glob
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.patches import Patch
import math
import numpy as np
from scipy.stats import t
from string import ascii_uppercase
from sklearn.metrics import r2_score
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression

TITLE_FONTSIZE = 18
AXIS_FONTSIZE = 14
FIG_WIDTH=12
FIG_HEIGHT=8
COLOR = {'lh_psp_n': 'C0', 'lh_psp_p': 'C0', 'rh_psp_n': 'C2', 'rh_psp_p': 'C2',
         'lh_fm_n': 'C1', 'lh_fm_p': 'C1', 'rh_fm_n': 'C3', 'rh_fm_p': 'C3'}

def setup_axis(ax, fontsize=5, ylim=[0, 0.04]):
    ax.set_xlim(0, 200)
    ax.set_ylim(ylim[0], ylim[1])
    ax.xaxis.set_major_locator(MultipleLocator(50))
    ax.yaxis.set_major_locator(MultipleLocator((ylim[1]-ylim[0])/4))
    for xv in [100, 150]:
        ax.axvline(x=xv, color='gray', linestyle='--', linewidth=0.8)
    ax.grid(True, which='major', linestyle='--', linewidth=0.5)
    ax.set_xlabel('Cortical depth (%)', fontsize=fontsize)
    ax.set_ylabel('Susceptibility (ppm)', fontsize=fontsize)

def plot_depthwise_profiles(hemi='rh', case=1):       
    ncols    = 10                     
    nrows    = math.ceil(n_subj / ncols) 
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols*1.5, nrows*1.5), sharex=True, sharey=True)
    if case == 1:
        prefix1, prefix2 = 'psp', 'fm'
    elif case == 2:
        prefix1, prefix2 = 'pc', 'fs'
    for ax, subj in zip(axes.flat, subjects):
        df= dfs[subj]
        df_std = dfs_std[subj]
        avg, std, ln = {}, {}, {}
        x  = df['depth'] * 10
        mask = (x >= 60)
        if hemi == 'lh' or hemi == 'rh':
            avg['psp_n'] = df[[f'{hemi}_{prefix1}_n']].mean(axis=1)
            avg['psp_p'] = df[[f'{hemi}_{prefix1}_p']].mean(axis=1)
            avg['fm_n']  = df[[f'{hemi}_{prefix2}_n']].mean(axis=1)
            avg['fm_p']  = df[[f'{hemi}_{prefix2}_p']].mean(axis=1)
            std['psp_n'] = df_std[[f'{hemi}_{prefix1}_n']].mean(axis=1)
            std['psp_p'] = df_std[[f'{hemi}_{prefix1}_p']].mean(axis=1)
            std['fm_n']  = df_std[[f'{hemi}_{prefix2}_n']].mean(axis=1)
            std['fm_p']  = df_std[[f'{hemi}_{prefix2}_p']].mean(axis=1)
        else:
            avg['psp_n'] = df[[f'lh_{prefix1}_n',f'rh_{prefix1}_n']].mean(axis=1)
            avg['psp_p'] = df[[f'lh_{prefix1}_p',f'rh_{prefix1}_p']].mean(axis=1)
            avg['fm_n']  = df[[f'lh_{prefix2}_n',f'rh_{prefix2}_n']].mean(axis=1)
            avg['fm_p']  = df[[f'lh_{prefix2}_p',f'rh_{prefix2}_p']].mean(axis=1)
            std['psp_n'] = df_std[[f'lh_{prefix1}_n',f'rh_{prefix1}_n']].mean(axis=1)
            std['psp_p'] = df_std[[f'lh_{prefix1}_p',f'rh_{prefix1}_p']].mean(axis=1)
            std['fm_n']  = df_std[[f'lh_{prefix2}_n',f'rh_{prefix2}_n']].mean(axis=1)
            std['fm_p']  = df_std[[f'lh_{prefix2}_p',f'rh_{prefix2}_p']].mean(axis=1)
        setup_axis(ax)
        ln['psp_n'], = ax.plot(x[mask], avg['psp_n'][mask], marker='o', label=f'{prefix1}_n', alpha=0.8, markersize=2)
        ln['psp_p'], = ax.plot(x[mask], avg['psp_p'][mask], marker='s', linestyle='--', label=f'{prefix1}_p', alpha=0.8, markersize=2)
        ln['fm_n'],  = ax.plot(x[mask], avg['fm_n'][mask], marker='^', label=f'{prefix2}_n', alpha=0.8, markersize=2)
        ln['fm_p'],  = ax.plot(x[mask], avg['fm_p'][mask], marker='v', linestyle='--', label=f'{prefix2}_p', alpha=0.8, markersize=2)
        # std shading
        for measure in ['psp_n', 'psp_p', 'fm_n', 'fm_p']:
            ax.fill_between(
                x[mask],
                (avg[measure] - std[measure])[mask],
                (avg[measure] + std[measure])[mask],
                color=ln[measure].get_color(), alpha=0.1
            )
        ax.set_title(subj, fontsize='x-small')
    for ax in axes.flat[n_subj:]:
        ax.set_visible(False)
    axes_list = axes.flatten()
    first_ax = axes_list[0]
    handles, labels = first_ax.get_legend_handles_labels()
    fig.legend(
        handles, labels,
        fontsize='xx-small', ncol=2,
        loc='upper right', bbox_to_anchor=(0.98,0.98)
    )
    fig.suptitle(f'All subjects — Depthwise Profiles {hemi} ({prefix1} & {prefix2})', fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()

def plot_avg_depthwise_profile(type='', case=1):
    if case == 1:
        prefix1, prefix2 = 'psp', 'fm'
        title1, title2 = 'precentral sulcus', 'middle frontal sulcus'
    elif case == 2:
        prefix1, prefix2 = 'pc', 'fs'
        title1, title2 = 'postcentral sulcus', 'superior frontal sulcus'
    x = dfs['Subj501']['depth']*10
    mask = (x >= 60)
    fig, (ax_psp, ax_fm) = plt.subplots(
        nrows=1, ncols=2,
        figsize=(FIG_WIDTH, FIG_HEIGHT),
        sharex=True,
        constrained_layout=True
    )
    setup_axis(ax_psp, fontsize=AXIS_FONTSIZE)
    setup_axis(ax_fm, fontsize=AXIS_FONTSIZE)
    list, df, avg, ln = {}, {}, {}, {}
    if type == 'hemi':
        valid_subj = {}
        for measure in [f'{prefix1}_n', f'{prefix1}_p', f'{prefix2}_n', f'{prefix2}_p']:
            valid_subj[measure] = [
                subj for subj, df in dfs.items()
                if df[f'lh_{measure}'].notnull().all() and df[f'rh_{measure}'].notnull().all()
            ]
            print(len(valid_subj[measure]))

        # ── sulcus 1 (prefix1) ────────────────────────────────────────────────
        for measure in [f'lh_{prefix1}_n', f'lh_{prefix1}_p', f'rh_{prefix1}_n', f'rh_{prefix1}_p']:
            # 값 평균(hemi별)
            key = f'{prefix1}_n' if measure.endswith('_n') else f'{prefix1}_p'
            series_list   = [ dfs[subj][measure]     for subj in valid_subj[key] ]
            df[measure]   = pd.concat(series_list, axis=1)             # rows=depth, cols=subjects
            avg[measure]  = df[measure].mean(axis=1)
            ln[measure], = ax_psp.plot(x[mask], avg[measure][mask], marker='o', label=measure, color=COLOR[measure])

            # PROP band: within-subject SDs 전파 (dfs_std에서 동일 measure 사용)
            sd_series_list = [ dfs_std[subj][measure] for subj in valid_subj[key] ]
            sd_df          = pd.concat(sd_series_list, axis=1)         # rows=depth, cols=subjects (s_i)
            n_subj         = len(valid_subj[key])
            band           = np.sqrt((sd_df ** 2).sum(axis=1)) / n_subj  # SD(mean) = sqrt(sum s_i^2)/n

            ax_psp.fill_between(
                x[mask],
                (avg[measure] - band)[mask],
                (avg[measure] + band)[mask],
                color=ln[measure].get_color(), alpha=0.15
            )

        # ── sulcus 2 (prefix2) ────────────────────────────────────────────────
        for measure in [f'lh_{prefix2}_n', f'lh_{prefix2}_p', f'rh_{prefix2}_n', f'rh_{prefix2}_p']:
            key = f'{prefix2}_n' if measure.endswith('_n') else f'{prefix2}_p'
            series_list   = [ dfs[subj][measure]     for subj in valid_subj[key] ]
            df[measure]   = pd.concat(series_list, axis=1)
            avg[measure]  = df[measure].mean(axis=1)
            ln[measure], = ax_fm.plot(x[mask], avg[measure][mask], marker='o', label=measure, color=COLOR[measure])

            sd_series_list = [ dfs_std[subj][measure] for subj in valid_subj[key] ]
            sd_df          = pd.concat(sd_series_list, axis=1)
            n_subj         = len(valid_subj[key])
            band           = np.sqrt((sd_df ** 2).sum(axis=1)) / n_subj

            ax_fm.fill_between(
                x[mask],
                (avg[measure] - band)[mask],
                (avg[measure] + band)[mask],
                color=ln[measure].get_color(), alpha=0.15
            )

        for i, ax in enumerate([ax_psp, ax_fm]):
            ax.text(0.01, 0.99, f'({ascii_uppercase[i]})', transform=ax.transAxes, ha='left', va='top', fontsize=AXIS_FONTSIZE)
            

    else:
        BAND = 'PROP'   # 'SEM' (권장) 또는 'PROP' 중 선택
        measures = [f'{prefix1}_n', f'{prefix1}_p', f'{prefix2}_n', f'{prefix2}_p']
        n = len(dfs)  # number of subjects

        vals_map   = {}   # per-measure: DataFrame (rows=depth, cols=subjects) of values
        mean_map   = {}   # per-measure: mean across subjects at each depth
        band_map   = {}   # per-measure: half-width of the shaded band to plot

        # (1) 값과 평균 계산
        for measure in measures:
            # subject별로 (lh, rh) 평균 후, 열 방향으로 concat → cols = subjects
            subj_series = [ df[[f'lh_{measure}', f'rh_{measure}']].mean(axis=1) for df in dfs.values() ]
            vals_df     = pd.concat(subj_series, axis=1)
            vals_map[measure] = vals_df
            mean_map[measure] = vals_df.mean(axis=1)

        # (2) 밴드(표준오차) 계산
        for measure in measures:
            if BAND == 'SEM':
                # across-subject SD (between-subject variability)
                sd_between = vals_map[measure].std(axis=1, ddof=1)
                band = sd_between / np.sqrt(n)  # SEM = s/sqrt(n)
            elif BAND == 'PROP':
                # 오차 전파: 각 subject의 within-subject SD(추정오차)가 dfs_std에 있음
                subj_sd_series = [
                    df_std[[f'lh_{measure}', f'rh_{measure}']].mean(axis=1)
                    for df_std in dfs_std.values()
                ]
                sd_df = pd.concat(subj_sd_series, axis=1)  # rows=depth, cols=subjects (y_i)
                band  = np.sqrt((sd_df ** 2).sum(axis=1)) / n  # SD(mean) = sqrt(sum y_i^2)/n
            else:
                raise ValueError("BAND must be 'SEM' or 'PROP'")
            band_map[measure] = band

        # (3) 플롯
        for measure in measures:
            ax = ax_psp if (measure.startswith(prefix1)) else ax_fm
            y  = mean_map[measure]
            b  = band_map[measure]
            ln[measure], = ax.plot(x[mask], y[mask], marker='o', label=measure)
            ax.fill_between(
                x[mask],
                (y - b)[mask],
                (y + b)[mask],
                color=ln[measure].get_color(), alpha=0.15
            )
    ax_psp.set_title(title1, fontsize=TITLE_FONTSIZE)
    ax_fm.set_title(title2, fontsize=TITLE_FONTSIZE)
    ax_psp.legend()
    ax_fm.legend()
    if type == 'hemi':
        plt.savefig('figures/avg_depthwise_profile_hemi.png')
    else:
        plt.savefig('figures/avg_depthwise_profile.png')

    # plt.tight_layout()
    plt.show()

def plot_age_depthwise_intensity(hemi, measure):
    depths = list(range(60, 210, 10))
    slopes, a0_list, a1_list, a2_list = [], [], [], []
    ncols = 5
    nrows = 3
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols*3, nrows*2.5), sharex=True, sharey=True)
    fig2, ax2 = plt.subplots(1, 1, figsize=(6, 8))
    # setup_axis(ax2, fontsize=12, ylim=[-0.0002, 0.0002])
    axes = axes.flatten()

    for ax, depth in zip(axes, depths):
        ax.set_ylim(0, 0.04)
        ages, vals, colors = [], [], []
        for subj in subjects:
            df = dfs[subj]
            age = df['Age'].iloc[0]
            x = df['depth'] * 10
            mask = (x == depth)

            if hemi == 'lh':
                col = f'lh_{measure}'
                val = df.loc[mask, col].iat[0]
            elif hemi == 'rh':
                col = f'rh_{measure}'
                val = df.loc[mask, col].iat[0]
            else:
                cols = [f'lh_{measure}', f'rh_{measure}']
                val = df.loc[mask, cols].mean(axis=1).iat[0]

            gender = df['Gender'].iloc[0]
            if gender.upper().startswith('F'):
                colors.append('red')
            elif gender.upper().startswith('M'):
                colors.append('blue')
            else:
                colors.append('black')
            ages.append(age)
            vals.append(val)

        ax.scatter(ages, vals, c=colors, alpha=0.7, label='data', s=10)
        ages_arr = np.array(ages)
        vals_arr = np.array(vals)
        mask_valid = ~np.isnan(vals_arr)
        ages_valid = ages_arr[mask_valid]
        vals_valid = vals_arr[mask_valid]

        # 1차 회귀선
        slope, intercept = np.polyfit(ages_valid, vals_valid, 1)
        slopes.append(slope)
        print(f"{slope:.3e}")
        ages_fit = np.linspace(min(ages), max(ages), 100)
        vals_fit = slope * ages_fit + intercept
        # 신뢰구간 계산
        n = len(ages_valid)
        vals_pred = slope * ages_valid + intercept
        r2 = r2_score(vals_valid, vals_pred)
        residuals = vals_valid - vals_pred
        sse = np.sum(residuals**2)
        std_err_est = np.sqrt(sse / (n - 2))
        t_critical = t.ppf((1 + 0.95) / 2, df=n - 2)
        margin_of_error = t_critical * std_err_est * np.sqrt(1/n + (ages_fit - np.mean(ages_valid))**2 / np.sum((ages_valid - np.mean(ages_valid))**2))
        ci_lower = vals_fit - margin_of_error
        ci_upper = vals_fit + margin_of_error
        ax.plot(ages_fit, vals_fit, label='regression', linewidth = 2)
        # ax.text(0.05, 0.95, f'$R^2 = {r2:.3f}$', transform=ax.transAxes, ha='left', va='top')
        ax.plot(ages_fit, ci_lower, color='black', linestyle='--', linewidth=1, label='95% Confidence Interval')
        ax.plot(ages_fit, ci_upper, color='black', linestyle='--', linewidth=1)
        # labels = [s for s in subjects]
        # for x_i, y_i, lab in zip(ages, vals, labels):
        #     ax.annotate(lab, (x_i, y_i), textcoords="offset points", xytext=(0, 4), ha='center', va='bottom', fontsize=5)
        ax.set_title(f'{depth}% depth', fontsize=12)
        ax.set_xlabel('Age (year)')
        ax.set_ylabel('Intensity (ppm)')

    for ax in axes[len(depths):]:
        ax.set_visible(False)
    legend_handles = [
        Patch(color='red',   label='Female'),
        Patch(color='blue',  label='Male'),
    ]
    fig.legend(handles=legend_handles, title='Gender', loc='upper right', fontsize=7.5)

    ax2.plot(depths, slopes, marker='o', linestyle='None', label='slopes')
    for xv in [100, 150]:
        ax2.axvline(x=xv, color='gray', linestyle='--', linewidth=0.8)
    ax2.axhline(0, color='gray', linestyle='--', linewidth=0.8)
    # ln_a2, = ax2.plot(depths, a2_list, marker='o', label='a2 (quadratic term)')
    # ln_a1, = ax2.plot(depths, a1_list, marker='s', label='a1 (linear term)')
    # ln_a0, = ax2.plot(depths, a0_list, marker='^', label='a0 (intercept)')

    ax2.set_xlabel('Depth (%)')
    ax2.set_ylabel('Coefficient value')
    ax2.set_title('slopes vs depth')
    ax2.legend()
    fig.suptitle(f'{measure} vs Age at each depth (hemi={hemi})', fontsize=TITLE_FONTSIZE)
    fig.tight_layout(rect=[0,0,1,0.96])
    fig.savefig('figures/agewise_intensity')
    plt.show()

def plot_age_average_intensity(hemi, measure):
    fig, ax = plt.subplots(1, 1, figsize=(9, 7.5))
    # fig2, ax2 = plt.subplots(1, 1, figsize=(9, 7.5))
    ax.set_ylim(0, 0.04)
    ages, vals, vals_std, colors, subs = [], [], [], [], []
    for subj in subjects:
        df = dfs[subj]
        age = df['Age'].iloc[0]
        x   = df['depth'] * 10
        mask = (x >= 60)
        if hemi == 'lh':
            col = f'lh_{measure}'
            val = df.loc[mask, col].mean()
        elif hemi == 'rh':
            col = f'rh_{measure}'
            val = df.loc[mask, col].mean()
        else:
            cols = [f'lh_{measure}', f'rh_{measure}']
            val = df.loc[mask, cols].mean(axis=1).mean()

        gender = df['Gender'].iloc[0]
        if gender.upper().startswith('F'):
            colors.append('red')
        elif gender.upper().startswith('M'):
            colors.append('blue')
        else:
            colors.append('black')
        ages.append(age)
        vals.append(val)
        subs.append(subj.replace('Subj', ''))
    ax.scatter(ages, vals, c=colors, alpha=0.7, label='data', s=10)
    ages_arr = np.array(ages)
    vals_arr = np.array(vals)
    mask_valid  = ~np.isnan(vals_arr)
    ages_valid  = ages_arr[mask_valid]
    vals_valid  = vals_arr[mask_valid]

    X = ages_valid.reshape(-1, 1)

    # Design matrix: [1, x, x^2]
    poly   = PolynomialFeatures(degree=2, include_bias=True)
    X_poly = poly.fit_transform(X)

    # Fit without separate intercept (bias column included)
    reg = LinearRegression(fit_intercept=False).fit(X_poly, vals_valid)
    a0, a1, a2 = reg.coef_[0], reg.coef_[1], reg.coef_[2]

    # 콘솔에 coeffs 출력 (원하면 depth 등 정보도 함께 출력)
    print(f"a2={a2:.3e}, a1={a1:.3e}, a0={a0:.3e}")

    # Smooth curve
    ages_fit = np.linspace(ages_valid.min(), ages_valid.max(), 200)
    X_fit    = poly.transform(ages_fit.reshape(-1, 1))
    vals_fit = X_fit @ reg.coef_

    # R^2
    y_pred_valid = X_poly @ reg.coef_
    r2 = r2_score(vals_valid, y_pred_valid)

    # 95% CI for mean prediction
    n = len(ages_valid)
    p = X_poly.shape[1]      # 3
    residuals = vals_valid - y_pred_valid
    sse   = np.sum(residuals**2)
    sigma2 = sse / (n - p)
    XtX_inv = np.linalg.inv(X_poly.T @ X_poly)
    se_fit = np.sqrt(np.sum((X_fit @ XtX_inv) * X_fit, axis=1) * sigma2)
    t_crit = t.ppf(0.975, df=n - p)
    ci_lower = vals_fit - t_crit * se_fit
    ci_upper = vals_fit + t_crit * se_fit

    # Plot on ax
    reg_line, = ax.plot(ages_fit, vals_fit, linewidth=2, label=f'Quadratic (R²={r2:.3f})')
    ax.plot(ages_fit, ci_lower, color='black', linestyle='--', linewidth=1)
    ax.plot(ages_fit, ci_upper, color='black', linestyle='--', linewidth=1)


    ax.set_xlabel('Age (year)')
    ax.set_ylabel('Average Intensity (ppm)')
    ax.set_title(f'{hemi} {measure} age vs average intensity')
    labels = [s for s in subs]
    for x_i, y_i, lab in zip(ages, vals, labels):
        ax.annotate(lab, (x_i, y_i), textcoords="offset points", xytext=(0, 4), ha='center', va='bottom', fontsize=7)
    legend_handles = [
        Patch(color='red',   label='Female'),
        Patch(color='blue',  label='Male'),
    ]
    ax.legend(handles=legend_handles, title='Gender')

    ax.text(0.98, 0.02, f'$R^2={r2:.3f}$', transform=ax.transAxes, ha='right', va='bottom', fontsize=9, bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=0.7))
    # for subj in subjects:
    #     df = dfs_std[subj]
    #     age = df['Age'].iloc[0]
    #     x   = df['depth'] * 10
    #     mask = (x >= 60)
    #     if hemi == 'lh':
    #         col = f'lh_{measure}'
    #         val = df.loc[mask, col].mean()
    #     elif hemi == 'rh':
    #         col = f'rh_{measure}'
    #         val = df.loc[mask, col].mean()
    #     else:
    #         cols = [f'lh_{measure}', f'rh_{measure}']
    #         val = df.loc[mask, cols].mean(axis=1).mean()
    #     vals_std.append(val)
    # ax2.scatter(ages, vals_std, c=colors, alpha=0.7, label='data', s=10)
    # # 회귀선
    # ages_arr = np.array(ages)
    # vals_std_arr = np.array(vals_std)
    # mask_valid = ~np.isnan(vals_std_arr)
    # ages_valid = ages_arr[mask_valid]
    # vals_std_valid = vals_std_arr[mask_valid]
    # slope, intercept = np.polyfit(ages_valid, vals_std_valid, 1)
    # ages_fit = np.linspace(min(ages), max(ages), 100)
    # vals_std_fit = slope * ages_fit + intercept
    # # 신뢰구간 계산
    # n = len(ages_valid)
    # vals_std_pred = slope * ages_valid + intercept
    # residuals = vals_std_valid - vals_std_pred
    # sse = np.sum(residuals**2)
    # std_err_est = np.sqrt(sse / (n - 2))
    # t_critical = t.ppf((1 + 0.95) / 2, df=n - 2)
    # margin_of_error = t_critical * std_err_est * np.sqrt(1/n + (ages_fit - np.mean(ages_valid))**2 / np.sum((ages_valid - np.mean(ages_valid))**2))
    # ci_lower = vals_std_fit - margin_of_error
    # ci_upper = vals_std_fit + margin_of_error
    # ax2.plot(ages_fit, vals_std_fit, label='regression', linewidth = 2)
    # ax2.plot(ages_fit, ci_lower, color='black', linestyle='--', linewidth=1, label='95% Confidence Interval')
    # ax2.plot(ages_fit, ci_upper, color='black', linestyle='--', linewidth=1)
    # for x_i, y_i, lab in zip(ages, vals_std, labels):
    #     ax2.annotate(lab, (x_i, y_i), textcoords="offset points", xytext=(0, 4), ha='center', va='bottom', fontsize=7)
    # ax2.set_xlabel('Age (year)')
    # ax2.set_ylabel('Average Standard Deviation (ppm)')
    # ax2.set_title(f'{hemi} {measure} age vs average std')
    # legend_handles = [
    #     Patch(color='red',   label='Female'),
    #     Patch(color='blue',  label='Male'),
    # ]
    # ax2.legend(handles=legend_handles, title='Gender')
    plt.show()

def plot_intensity_roundage_rounddepth(type=['mean']):
    # --- normalize `type` -> list[str] ---
    if type is None:
        types = ['all']
    elif isinstance(type, str):
        types = [type.lower()]
    else:
        # iterable of strings
        types = [str(t).lower() for t in type]

    df_irard = pd.read_excel(os.path.join('intensity_roundage_rounddepth.xlsx'))
    x = dfs['Subj035']['depth'] * 10
    mask = (x >= 60)

    fig, (ax_n, ax_p) = plt.subplots(
        nrows=1, ncols=2,
        figsize=(FIG_WIDTH, FIG_HEIGHT),
        sharex=True,
        constrained_layout=True
    )

    setup_axis(ax_n, fontsize=AXIS_FONTSIZE, ylim=[-0.0002, 0.0002])
    setup_axis(ax_p, fontsize=AXIS_FONTSIZE, ylim=[-0.00008, 0.00008])
    for ax in (ax_n, ax_p):
        ax.set_axisbelow(True)
        ax.axhline(0, color='0.2', linewidth=1.6, alpha=0.9, zorder=0)

    ax_n.set_ylabel('Susceptibility / Age (ppm/year)', fontsize=AXIS_FONTSIZE)
    ax_p.set_ylabel('Susceptibility / Age (ppm/year)', fontsize=AXIS_FONTSIZE)
    ax_p.ticklabel_format(axis='y', style='plain')

    # --- group definitions ---
    groups = {
        'psp': ['lh_psp_n', 'lh_psp_p', 'rh_psp_n', 'rh_psp_p'],
        'fm' : ['lh_fm_n',  'lh_fm_p',  'rh_fm_n',  'rh_fm_p' ],
        'pc' : ['lh_pc_n',  'lh_pc_p',  'rh_pc_n',  'rh_pc_p' ],
        'fs' : ['lh_fs_n',  'lh_fs_p',  'rh_fs_n',  'rh_fs_p' ],
    }
    all_measures = sum(groups.values(), [])  # flatten

    did_any = False

    # 1) mean block (independent)
    if 'mean' in types:
        n_measures = [m for m in all_measures if m.endswith('_n')]
        p_measures = [m for m in all_measures if m.endswith('_p')]

        mean_n = df_irard[n_measures].mean(axis=1)
        mean_p = df_irard[p_measures].mean(axis=1)
        ax_n.plot(x[mask], mean_n[mask], marker='o', label='n measures Mean')
        ax_p.plot(x[mask], mean_p[mask], marker='o', label='p measures Mean')
        did_any = True

    # 2) requested ROI groups (can be multiple)
    requested_groups = {g for g in groups.keys() if g in types}
    if requested_groups:
        for g in requested_groups:
            for measure in groups[g]:
                if measure.endswith('_n'):
                    ax_n.plot(x[mask], df_irard[measure][mask], marker='o',
                              label=measure, linestyle='--', color=COLOR[measure])
                else:
                    ax_p.plot(x[mask], df_irard[measure][mask], marker='o',
                              label=measure, linestyle='--', color=COLOR[measure])
        did_any = True

    # 3) default or 'all' → plot everything
    if ('all' in types) or (not did_any):
        for measure in all_measures:
            if measure.endswith('_n'):
                ax_n.plot(x[mask], df_irard[measure][mask], marker='o',
                          label=measure, linestyle='none')
            else:
                ax_p.plot(x[mask], df_irard[measure][mask], marker='o',
                          label=measure, linestyle='none')

    ax_n.set_title('Depth profile of age-related changes \n in negative susceptibiliy', fontsize=AXIS_FONTSIZE)
    ax_p.set_title('Depth profile of age-related changes \n in positive susceptibility', fontsize=AXIS_FONTSIZE)
    ax_n.legend(loc='lower left')
    ax_p.legend(loc='lower left')
    plt.savefig('figures/intensity_roundage_rounddepth.png', bbox_inches='tight')
    plt.show()

def fit_irard_slope_xint(measure: str,
                         depth_start: float,
                         depth_end: float,
                         excel_path: str = 'intensity_roundage_rounddepth.xlsx',
                         subject_for_depth: str | None = None):
    df_irard = pd.read_excel(os.path.join(excel_path))

    if subject_for_depth is None:
        if len(dfs) == 0:
            raise RuntimeError("dfs is empty — cannot build depth vector.")
        subject_for_depth = next(iter(dfs.keys()))
    depth_vec = (dfs[subject_for_depth]['depth'] * 10).to_numpy()  # % depth

    # 3) 유효 구간 마스크
    mask = (depth_vec >= depth_start) & (depth_vec <= depth_end)
    if measure not in df_irard.columns:
        raise KeyError(f"measure '{measure}' not in df_irard. Available: {list(df_irard.columns)}")

    y_all = df_irard[measure].to_numpy()
    x_seg = depth_vec[mask]
    y_seg = y_all[mask]

    # 결측 제거
    valid = ~np.isnan(y_seg)
    x = x_seg[valid].reshape(-1, 1)
    y = y_seg[valid]

    # 4) 선형 회귀 (y = m*x + b)
    reg = LinearRegression(fit_intercept=True).fit(x, y)
    slope = float(reg.coef_[0])
    intercept = float(reg.intercept_)

    # 5) x-절편 (y=0 → x = -b/m)
    x_intercept = (-intercept / slope) if not np.isclose(slope, 0.0) else np.nan

    print(f"[{measure}] depth {depth_start:.0f}–{depth_end:.0f}% "
          f"→ slope = {slope:.3e} (ppm/year per %depth), x-intercept = {x_intercept:.2f}%")
    return slope, x_intercept

def plot_bar_slope_xint(xlsx_path='slope_xint.xlsx', savepath=None, rotate=45):
    """
    2x2 bar charts:
      [0,0] slope (suffix '_n'), [0,1] slope ('_p')
      [1,0] intercept ('_n'),    [1,1] intercept ('_p')

    Colors: (psp + pc) same, (fm + fs) same.
    """
    # 1) load & normalize labels
    df = pd.read_excel(xlsx_path)
    first_col = df.columns[0]
    df.set_index(first_col, inplace=True)
    df.index   = df.index.map(lambda s: str(s).strip().lower())
    df.columns = [str(c).strip().lower() for c in df.columns]

    # 행 라벨 찾기(둘 다 지원: 'intercept' 또는 'intercept (%)')
    def find_row(target_kw: str) -> str:
        if target_kw in df.index:
            return target_kw
        cand = [idx for idx in df.index if target_kw in idx]
        if cand:
            return cand[0]
        raise ValueError(f"Row containing '{target_kw}' not found. Rows: {list(df.index)}")

    row_slope     = find_row('slope')
    row_intercept = find_row('intercept')

    # 2) helpers
    order_map   = {'psp':0, 'pc':1, 'fm':2, 'fs':3}
    hemi_map    = {'lh':0, 'rh':1}
    colors_map  = {'psp':'tab:blue', 'pc':'tab:blue', 'fm':'tab:orange', 'fs':'tab:orange'}

    def parse_tokens(col):
        # ex) 'lh_psp_n' → hemi='lh', area='psp', suf='n'
        parts = col.split('_')
        if len(parts) != 3:
            return None, None, None
        hemi, area, suf = parts
        return hemi, area, suf

    def ordered_cols(suffix):  # suffix: 'n' or 'p'
        cols = [c for c in df.columns if c.endswith('_'+suffix)]
        cols = [c for c in cols if parse_tokens(c)[1] in order_map and parse_tokens(c)[0] in hemi_map]
        cols.sort(key=lambda c: (order_map[parse_tokens(c)[1]], hemi_map[parse_tokens(c)[0]]))
        return cols

    def bar_on_ax(ax, rowname, suffix, title):
        cols = ordered_cols(suffix)
        if len(cols) == 0:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center')
            return

        vals = df.loc[rowname, cols].astype(float).to_numpy()
        x = np.arange(len(cols))

        # NaN은 0으로 두고, 이후 완전 투명 처리
        heights = np.nan_to_num(vals, nan=0.0)
        bars = ax.bar(x, heights, width=0.8, edgecolor='black', linewidth=0.6)

        for i, (bar, v) in enumerate(zip(bars, vals)):
            _, area, _ = parse_tokens(cols[i])
            if np.isnan(v):
                bar.set_facecolor((0, 0, 0, 0))
                bar.set_edgecolor((0, 0, 0, 0))
            else:
                bar.set_facecolor(colors_map[area])

        # 축 라벨(단위)
        if 'slope' in rowname:
            ylab = 'Slope (ppm/age per %depth)'
        else:
            ylab = 'Intercept (%)'

        ax.set_title(title, fontsize=TITLE_FONTSIZE)
        ax.set_xlabel(f"measure")
        ax.set_ylabel(ylab)
        ax.grid(True, axis='y', linestyle='--', linewidth=0.5, alpha=0.6)
        ax.axhline(0, color='0.3', lw=1.0, ls='--')
        labels = [c[: -len('_'+suffix)] for c in cols]  # 예: 'lh_psp_n' -> 'lh_psp'

        # usetex이면 '_'를 이스케이프
        if plt.rcParams.get('text.usetex', False):
            labels = [lab.replace('_', r'\_') for lab in labels]

        ax.set_xticks(np.arange(len(cols)))
        ax.set_xticklabels(labels, rotation=rotate, ha='right')

    # 3) figure
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)

    bar_on_ax(axes[0,0], row_slope,     'n', 'Slope (n)')
    bar_on_ax(axes[0,1], row_slope,     'p', 'Slope (p)')
    bar_on_ax(axes[1,0], row_intercept, 'n', 'Intercept (n)')
    bar_on_ax(axes[1,1], row_intercept, 'p', 'Intercept (p)')

    # legend: colors only by family
    handles = [Patch(color='tab:blue',   label='psp + pc'),
               Patch(color='tab:orange', label='fm + fs')]
    fig.legend(handles=handles, loc='upper right', frameon=False)

    fig.suptitle('Bar charts of slope & x-intercept', fontsize=TITLE_FONTSIZE)
    if savepath:
        fig.savefig(savepath, bbox_inches='tight', dpi=300)
    plt.show()


age_df = pd.read_csv( 
    'patient_ages.txt',
    sep=r'\s+',
    header=None,
    names=['SubjCode', 'Age', 'Gender', 'Extra'],
    usecols=[0,1,2]
)

excel_paths = glob.glob(os.path.join('Stat_final2', '*.xlsx'))
dfs = {
    os.path.splitext(os.path.basename(path))[0]: pd.read_excel(path, sheet_name='Avg')
    for path in excel_paths
}

dfs_std = {
    os.path.splitext(os.path.basename(path))[0]: pd.read_excel(path, sheet_name='Std')
    for path in excel_paths
}

excel_keys = set(dfs.keys())

# 기본 base: 'Subj' + 숫자 3자리. 엑셀 키에 'base_'가 있으면 '_'를 붙이고, 아니면 base 그대로.
age_df['Subject'] = (
    'Subj' + age_df['SubjCode'].str[1:].str.zfill(3)
).apply(lambda base: base + '_' if (base + '_') in excel_keys else base)

age_df = age_df.set_index('Subject')

for subj, df in dfs.items():
    df = df.rename(columns={'index': 'depth'})
    df['Age']    = age_df.loc[subj, 'Age']
    df['Gender'] = age_df.loc[subj, 'Gender']
    dfs[subj] = df
for subj, df in dfs_std.items():
    df = df.rename(columns={'index': 'depth'})
    df['Age']    = age_df.loc[subj, 'Age']
    df['Gender'] = age_df.loc[subj, 'Gender']
    dfs_std[subj] = df

subjects = list(dfs.keys())           
n_subj   = len(subjects)
print(n_subj)
# plot_depthwise_profiles(hemi='lh', case=2)
# plot_avg_depthwise_profile('hemi', case=1)
plot_age_depthwise_intensity(hemi='rh', measure='pc_n')
# plot_age_average_intensity(hemi='', measure='psp_n')
# plot_intensity_roundage_rounddepth(['psp', 'fm'])

# fit_irard_slope_xint('lh_psp_n', 100, 150)
# fit_irard_slope_xint('rh_psp_n', 100, 150)
# fit_irard_slope_xint('lh_fm_n', 100, 150)
# fit_irard_slope_xint('rh_fm_n', 100, 150)
# fit_irard_slope_xint('lh_pc_n', 100, 150)
# fit_irard_slope_xint('rh_pc_n', 100, 150)
# fit_irard_slope_xint('lh_fs_n', 100, 150)
# fit_irard_slope_xint('rh_fs_n', 100, 150)
# print()

# fit_irard_slope_xint('lh_psp_p', 110, 150)
# fit_irard_slope_xint('rh_psp_p', 100, 150)
# fit_irard_slope_xint('lh_fm_p', 120, 150)
# fit_irard_slope_xint('rh_fm_p', 120, 150)
# fit_irard_slope_xint('lh_pc_p', 120, 150)
# fit_irard_slope_xint('rh_pc_p', 130, 150)
# fit_irard_slope_xint('lh_fs_p', 110, 150)
# fit_irard_slope_xint('rh_fs_p', 110, 150)

# plot_bar_slope_xint('slope_xint.xlsx', savepath='figures/slope_xint_bars.png')
