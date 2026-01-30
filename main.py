import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import solve_ivp

# =========================
# tqdm（進捗バー）：無い環境でも落ちないようにする
# =========================
try:
    from tqdm.auto import tqdm
except Exception:
    tqdm = None

# =========================
# 定数 / 単位（SI）
# =========================
G  = 6.67430e-11
M_sun = 1.98847e30
MU = G * M_sun
C  = 299792458.0

DAY = 86400.0
YEAR = 365.25 * DAY
CENTURY = 100.0 * YEAR
RAD2ARCSEC = (180.0/np.pi) * 3600.0

# =========================
# 水星の代表軌道要素（概略）
# =========================
a = 57_909_050_000.0      # m
e = 0.205630

def state_at_perihelion(a, e, mu=MU):
    """近日点通過時の状態（簡易：軌道面をx-y平面、近日点を+X方向に配置）"""
    rp = a * (1.0 - e)
    vp = np.sqrt(mu * (1.0 + e) / (a * (1.0 - e)))
    r0 = np.array([rp, 0.0, 0.0])
    v0 = np.array([0.0, vp, 0.0])
    return r0, v0

# =========================
# 加速度モデル
# =========================
def a_newton(r, mu):
    """ニュートン重力（中心力）"""
    rn = np.linalg.norm(r)
    return -mu * r / rn**3

def a_1pn_harmonic(r, v, mu, c_eff):
    """
    1PN（調和座標）テスト粒子の標準形：
      a_1PN = (mu/(c^2 r^3)) * [ (4mu/r - v^2) r + 4 (r·v) v ]
    ※水星で ~43 arcsec/century を再現する用途
    """
    rn = np.linalg.norm(r)
    v2 = np.dot(v, v)
    rv = np.dot(r, v)
    return (mu / (c_eff**2 * rn**3)) * ((4.0*mu/rn - v2) * r + 4.0 * rv * v)

def a_2pn_eq1(r, v, mu, c_eff):
    """
    Iorio(2020) の Eq.(1)（direct 2PN 加速度、調和座標）：
      A_2PN = (mu^2)/(c^4 r^3) * [ (2 v_r^2 - 9 mu/r) r_hat - 2 v_r v ]
    ここで v_r = v · r_hat
    """
    rn = np.linalg.norm(r)
    rhat = r / rn
    vr = np.dot(v, rhat)
    return (mu**2 / (c_eff**4 * rn**3)) * ((2.0*vr**2 - 9.0*mu/rn) * rhat - 2.0*vr * v)

# =========================
# モデル切替（ワンタッチ）
# =========================
def make_rhs(model: str, mu: float, c_eff: float):
    """model: 'newton', 'pn1', 'pn2', 'pn12' を想定"""
    def rhs(t, y):
        r = y[:3]
        v = y[3:]
        acc = a_newton(r, mu)

        if model in ("pn1", "pn12"):
            acc = acc + a_1pn_harmonic(r, v, mu, c_eff)

        if model in ("pn2", "pn12"):
            acc = acc + a_2pn_eq1(r, v, mu, c_eff)

        return np.hstack([v, acc])
    return rhs

# =========================
# 近日点イベント：dr/dt = 0（近地点ではなく近日点側のみ拾う）
# =========================
def perihelion_event(t, y):
    r = y[:3]
    v = y[3:]
    rn = np.linalg.norm(r)
    return np.dot(r, v) / rn  # dr/dt

perihelion_event.terminal = False
perihelion_event.direction = +1

# =========================
# 近日点歳差率の推定（線形フィット）
# =========================
def fit_precession_arcsec_per_century(t_events, y_events):
    """
    近日点位置ベクトルの偏角 theta = atan2(y, x) を unwrap し、
    theta vs time(century) を一次近似して slope を得る
    """
    r = y_events[:, :3]
    theta = np.unwrap(np.arctan2(r[:, 1], r[:, 0]))

    tc = (t_events - t_events[0]) / CENTURY
    p = np.polyfit(tc, theta, 1)
    slope_rad_per_century = p[0]
    return slope_rad_per_century * RAD2ARCSEC, theta, tc

def theory_1pn_arcsec_per_century(a, e, mu, c_eff):
    """Einstein 1PN: Δω = 6πμ/(a(1-e^2)c^2) [rad/orbit] を century あたりに換算"""
    d_per_orbit = 6.0*np.pi*mu/(a*(1.0-e**2)*c_eff**2)
    P = 2.0*np.pi*np.sqrt(a**3/mu)
    orbits_per_century = CENTURY / P
    return d_per_orbit * orbits_per_century * RAD2ARCSEC

def theory_2pn_direct_arcsec_per_century(a, e, mu, c_eff):
    """
    direct 2PN の理論式（Iorio の direct 2PN、調和座標）
      ωdot_dir^(2PN) = n_b * μ^2 * (28 - e^2) / (4 c^4 a^2 (1-e^2)^2)
    を arcsec/century に換算
    """
    n_b = np.sqrt(mu / a**3)  # [rad/s]
    wdot = n_b * (mu**2) * (28.0 - e**2) / (4.0 * (c_eff**4) * (a**2) * (1.0 - e**2)**2)  # [rad/s]
    return wdot * CENTURY * RAD2ARCSEC

# =========================
# 根本対策：同一時刻で「近点方向」を評価して差分を取る
# =========================
def ecc_angle_on_times(sol, t_query, mu=MU):
    """
    共通時刻 t_query における偏心（LRL）ベクトルの向き（角度）を返す。
    これを「近点方向（オスキュレーティング近似）」として用いる。

    注意：
    - sol.sol(t) を使うので、solve_ivp は dense_output=True が前提
    """
    y = sol.sol(t_query)          # shape: (6, N)
    r = y[:3, :].T                # (N, 3)
    v = y[3:, :].T                # (N, 3)

    rn = np.linalg.norm(r, axis=1)[:, None]
    h  = np.cross(r, v)           # (N, 3)
    evec = np.cross(v, h) / mu - r / rn

    # 軌道面はほぼxy平面前提なので、x-yで角度を取る
    theta = np.unwrap(np.arctan2(evec[:, 1], evec[:, 0]))
    return theta

def diff_series_ecc(results, model_a, model_b, t_grid):
    """
    共通時刻 t_grid 上で、偏心ベクトル角の差分系列を作る。
    - 戻り値：x_year, dtheta_arcsec
    """
    th_a = ecc_angle_on_times(results[model_a]["sol"], t_grid, mu=MU)
    th_b = ecc_angle_on_times(results[model_b]["sol"], t_grid, mu=MU)
    d = (th_a - th_b)
    d = d - d[0]  # 初期オフセットを除去
    x_year = (t_grid - t_grid[0]) / YEAR
    return x_year, d * RAD2ARCSEC

# =========================
# 結果出力：CSV（seriesのみ）
# =========================
def export_results(prefix, summary_rows, df_series_dict):
    """
    結果をCSVに出力する（seriesのみ）
    - summaryはCSVにせず、別途 main() 側で .txt に出す
    - xlsxは出力しない
    """
    df_series = pd.DataFrame(df_series_dict)

    # CSV（Excel互換を意識してBOM付きUTF-8）
    df_series.to_csv(f"{prefix}_series.csv", index=False, encoding="utf-8-sig")

    print(f"[export] wrote: {prefix}_series.csv")
    return


# =========================
# 実行
# =========================
def run(model, years_long=100.0, scale_c=1.0, rtol=1e-12, atol=1e-15, max_step=20000.0):
    c_eff = C / scale_c
    r0, v0 = state_at_perihelion(a, e, MU)
    y0 = np.hstack([r0, v0])

    t_span = (0.0, years_long * YEAR)

    sol = solve_ivp(
        make_rhs(model, MU, c_eff),
        t_span,
        y0,
        method="DOP853",
        rtol=rtol,
        atol=atol,
        events=perihelion_event,
        max_step=max_step,
        dense_output=True   # 根本対策（同一時刻評価）に必要：既存の結果解釈は変えない
    )

    t_ev = sol.t_events[0]
    y_ev = sol.y_events[0]
    return sol, t_ev, y_ev, c_eff

def main():
    YEARS_LONG = 100.0
    SCALE_C = 1.0  # デバッグでPN効果を増幅したいなら増やす（例：100）

    models = ["newton", "pn1", "pn2", "pn12"]

    # 進捗バー（tqdmが無ければ普通のlist）
    model_iter = tqdm(models, desc="モデル計算", unit="model") if tqdm is not None else models

    results = {}
    for m in model_iter:
        sol, t_ev, y_ev, c_eff = run(m, years_long=YEARS_LONG, scale_c=SCALE_C)
        rate, theta, tc = fit_precession_arcsec_per_century(t_ev, y_ev)
        results[m] = dict(sol=sol, t_ev=t_ev, y_ev=y_ev, rate=rate, theta=theta, tc=tc, c_eff=c_eff)

    # Newton 基準差分
    base = results["newton"]["rate"]
    r_newton = results["newton"]["rate"] - base
    r_pn1    = results["pn1"]["rate"]    - base
    r_pn2    = results["pn2"]["rate"]    - base
    r_pn12   = results["pn12"]["rate"]   - base

    th_1pn = theory_1pn_arcsec_per_century(a, e, MU, results["pn1"]["c_eff"])
    th_2pn_dir = theory_2pn_direct_arcsec_per_century(a, e, MU, results["pn2"]["c_eff"])

    print("============================================================")
    print(f"SCALE_C         : {SCALE_C}  (c_eff = c / SCALE_C)")
    print(f"YEARS_LONG      : {YEARS_LONG}")
    print(f"perihelia used  : {len(results['newton']['t_ev'])}  (Newton)")
    print("------------------------------------------------------------")
    print("MODEL           : newton")
    print(f"fit rate        : {r_newton:.6f} arcsec/century  (difference vs Newton)")
    print("------------------------------------------------------------")
    print("MODEL           : pn1")
    print(f"fit rate        : {r_pn1:.6f} arcsec/century  (difference vs Newton)")
    print(f"theory 1PN      : {th_1pn:.6f} arcsec/century  (Einstein 1PN)")
    print(f"pn1 - theory    : {(r_pn1 - th_1pn)*1e6:.3f} microarcsec/century  (indirect-2PN scale check)")
    print("------------------------------------------------------------")
    print("MODEL           : pn2 (Eq.(1) only)")
    print(f"fit rate        : {r_pn2:.9f} arcsec/century  (difference vs Newton)")
    print(f"theory 2PN(dir) : {th_2pn_dir:.9f} arcsec/century  (direct 2PN theory)")
    print(f"pn2 - theory    : {(r_pn2 - th_2pn_dir)*1e6:.3f} microarcsec/century")
    print("------------------------------------------------------------")
    print("MODEL           : pn12 (1PN + Eq.(1))")
    print(f"fit rate        : {r_pn12:.6f} arcsec/century  (difference vs Newton)")
    print(f"pn12 - pn1      : {(r_pn12 - r_pn1)*1e6:.3f} microarcsec/century  (direct-2PN scale check)")
    print("============================================================")

    # =========================
    # 追加：コマンドプロンプト出力（サマリ）を .txt に保存
    # =========================
    EXPORT_PREFIX = "mercury_pn_validation"
    summary_txt_path = f"{EXPORT_PREFIX}_summary.txt"

    lines = []
    lines.append("============================================================")
    lines.append(f"SCALE_C         : {SCALE_C}  (c_eff = c / SCALE_C)")
    lines.append(f"YEARS_LONG      : {YEARS_LONG}")
    lines.append(f"perihelia used  : {len(results['newton']['t_ev'])}  (Newton)")
    lines.append("------------------------------------------------------------")
    lines.append("MODEL           : newton")
    lines.append(f"fit rate        : {r_newton:.6f} arcsec/century  (difference vs Newton)")
    lines.append("------------------------------------------------------------")
    lines.append("MODEL           : pn1")
    lines.append(f"fit rate        : {r_pn1:.6f} arcsec/century  (difference vs Newton)")
    lines.append(f"theory 1PN      : {th_1pn:.6f} arcsec/century  (Einstein 1PN)")
    lines.append(f"pn1 - theory    : {(r_pn1 - th_1pn)*1e6:.3f} microarcsec/century  (indirect-2PN scale check)")
    lines.append("------------------------------------------------------------")
    lines.append("MODEL           : pn2 (Eq.(1) only)")
    lines.append(f"fit rate        : {r_pn2:.9f} arcsec/century  (difference vs Newton)")
    lines.append(f"theory 2PN(dir) : {th_2pn_dir:.9f} arcsec/century  (direct 2PN theory)")
    lines.append(f"pn2 - theory    : {(r_pn2 - th_2pn_dir)*1e6:.3f} microarcsec/century")
    lines.append("------------------------------------------------------------")
    lines.append("MODEL           : pn12 (1PN + Eq.(1))")
    lines.append(f"fit rate        : {r_pn12:.6f} arcsec/century  (difference vs Newton)")
    lines.append(f"pn12 - pn1      : {(r_pn12 - r_pn1)*1e6:.3f} microarcsec/century  (direct-2PN scale check)")
    lines.append("============================================================")

    with open(summary_txt_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    print(f"[export] wrote: {summary_txt_path}")

    # =========================
    # 既存：xy 軌道プロット
    # =========================
    fig, ax = plt.subplots(figsize=(8, 6))
    for m in ["newton", "pn1", "pn2", "pn12"]:
        sol = results[m]["sol"]
        x = sol.y[0]
        y = sol.y[1]
        ax.plot(x/1e9, y/1e9, label=m)

    ax.set_aspect("equal", "box")
    ax.set_xlabel("x [Gm]")
    ax.set_ylabel("y [Gm]")
    ax.grid(True)

    # 凡例を外側へ
    ax.legend(bbox_to_anchor=(1.02, 1.0), loc="upper left", borderaxespad=0.0)
    plt.tight_layout(rect=[0.0, 0.0, 0.8, 1.0])
    plt.show()

    # =========================
    # 根本対策版：共通時刻で差分を作る
    # - 近日点イベント時刻列（Newton側）を共通グリッドとして採用
    #   （1周に1点なので、長期トレンド確認にちょうど良い）
    # =========================
    t_grid = results["newton"]["t_ev"]

    # =========================
    # 図2：近日点方向角の差（pn1 - newton）[arcsec vs year]
    # =========================
    x2, y2 = diff_series_ecc(results, "pn1", "newton", t_grid)

    fig2, ax2 = plt.subplots(figsize=(8, 5))
    ax2.plot(x2, y2, label="pn1 - newton")
    ax2.set_xlabel("year")
    ax2.set_ylabel("perihelion direction change [arcsec]")
    ax2.grid(True)
    ax2.legend(bbox_to_anchor=(1.02, 1.0), loc="upper left", borderaxespad=0.0)
    plt.tight_layout(rect=[0.0, 0.0, 0.8, 1.0])
    plt.show()

    # =========================
    # 図3：direct 2PNの見比べ（pn2-newton と pn12-pn1）[μas vs year]
    # =========================
    x3a, y3a = diff_series_ecc(results, "pn2",  "newton", t_grid)
    x3b, y3b = diff_series_ecc(results, "pn12", "pn1",    t_grid)

    fig3, ax3 = plt.subplots(figsize=(8, 5))
    ax3.plot(x3a, y3a * 1e6, label="pn2 - newton")
    ax3.plot(x3b, y3b * 1e6, label="pn12 - pn1")
    ax3.set_xlabel("year")
    ax3.set_ylabel("perihelion direction change [microarcsec]")
    ax3.grid(True)
    ax3.legend(bbox_to_anchor=(1.02, 1.0), loc="upper left", borderaxespad=0.0)
    plt.tight_layout(rect=[0.0, 0.0, 0.8, 1.0])
    plt.show()

    # =========================
    # 追加：series.csv 出力（summaryはtxtで出しているので不要）
    # =========================
    EXPORT_PREFIX = "mercury_pn_validation"

    df_series_dict = {
        "year": x3a,
        "pn1_minus_newton_arcsec": y2,
        "pn2_minus_newton_microarcsec": y3a * 1e6,
        "pn12_minus_pn1_microarcsec": y3b * 1e6,
    }

    # summary_rows は export_results では使わないが、引数互換のためダミーで渡す
    export_results(EXPORT_PREFIX, summary_rows=[], df_series_dict=df_series_dict)



if __name__ == "__main__":
    main()
