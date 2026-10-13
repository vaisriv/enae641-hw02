from functools import reduce
from pathlib import Path

import control as ct
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad, quad_vec, solve_ivp
from scipy.linalg import expm

matplotlib.use("Agg")

ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "outputs/figures"
TEXT = ROOT / "outputs/text"

# series colours (first three categorical slots, validated as a set) and chart ink
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED = "#0b0b0b", "#52514e"
plt.rcParams.update(
    {
        "axes.edgecolor": "#c3c2b7",
        "axes.labelcolor": INK,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.titlesize": 9,
        "font.size": 9,
        "grid.color": "#e1e0d9",
        "grid.linewidth": 0.6,
        "legend.frameon": False,
        "lines.linewidth": 1.6,
        "text.color": INK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
    }
)


def heading(problem):
    print(f"\n------\n {problem}\n------")


def report(name, content):
    print(content)
    (TEXT / f"{name}.txt").write_text(content + "\n", encoding="utf-8")


def save(fig, name):
    fig.tight_layout()
    fig.savefig(FIGURES / f"{name}.png", dpi=200)
    plt.close(fig)


def show(M, precision=4):
    """Format an array, rounding away floating-point dust (and negative zeros)."""
    return np.array2string(
        np.round(np.asarray(M), 10) + 0.0, precision=precision, suppress_small=True
    )


def poly_str(coeffs, var="s"):
    """Format polynomial coefficients (highest power first) as text."""
    coeffs = np.trim_zeros(np.round(np.real(coeffs), 6) + 0.0, "f")
    terms = []
    for power, c in zip(range(coeffs.size - 1, -1, -1), coeffs):
        if c == 0:
            continue
        mag = "" if abs(c) == 1 and power > 0 else f"{abs(c):g}"
        sym = "" if power == 0 else var if power == 1 else f"{var}^{power}"
        terms.append(("-" if c < 0 else "+", f"{mag} {sym}".strip()))
    if not terms:
        return "0"
    text = ("-" if terms[0][0] == "-" else "") + terms[0][1]
    return text + "".join(f" {sign} {body}" for sign, body in terms[1:])


def tf_str(g):
    """Format a SISO transfer function with a monic denominator."""
    num, den = np.asarray(g.num[0][0]), np.asarray(g.den[0][0])
    num, den = num / den[0], den / den[0]
    numerator = poly_str(num)
    if numerator == "0":
        return "0"
    if any(op in numerator[1:] for op in "+-"):
        numerator = f"({numerator})"
    return f"{numerator} / ({poly_str(den)})"


def channel_tfs(A, B, C, cancel=True):
    """Transfer function of every input-output channel of the model (A, B, C, 0)."""
    G = {}
    for i in range(C.shape[0]):
        for j in range(B.shape[1]):
            g = ct.tf(ct.ss(A, B[:, [j]], C[[i], :], 0))
            num = np.array(g.num[0][0], dtype=float)
            num[np.abs(num) < 1e-9] = 0.0  # round-off in an identically zero channel
            g = ct.tf(num, g.den[0][0]) if num.any() else ct.tf([0.0], [1.0])
            G[i, j] = g.minreal(tol=1e-6) if cancel else g
    return G


def stm_numeric(A, tau, t):
    """Integrate dPhi/dt = A(t) Phi from Phi(tau, tau) = I, returning Phi(t_k, tau)."""
    n = A(tau).shape[0]
    sol = solve_ivp(
        lambda s, phi: (A(s) @ phi.reshape(n, n)).ravel(),
        (tau, t[-1]),
        np.eye(n).ravel(),
        method="DOP853",
        t_eval=t,
        rtol=1e-11,
        atol=1e-12,
    )
    assert sol.success
    return sol.y.T.reshape(-1, n, n)


def stm_error(A, Phi, taus=(-1.0, 0.0, 1.5), span=3.0):
    """Largest relative mismatch between a closed-form STM and the numerical STM."""
    worst = 0.0
    for tau in taus:
        for end in (tau + span, tau - span):
            t = np.linspace(tau, end, 31)
            numeric = stm_numeric(A, tau, t)
            exact = np.array([Phi(tk, tau) for tk in t])
            error = np.abs(numeric - exact) / np.maximum(1.0, np.abs(exact))
            worst = max(worst, error.max())
    assert worst < 1e-7
    return worst


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)
    TEXT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(641)

    #########
    # p01 #
    #########
    heading("p01")
    A = np.array([[0, 0, -13], [1, 0, -17], [0, 1, -5]], dtype=float)
    B = np.array([[13, 1], [4, 1], [1, 0]], dtype=float)
    C = np.array([[1, -1, 1], [0, 1, -1]], dtype=float)

    # a) eigenvectors from the companion structure, v(lambda) = [l^2 + 5l + 17, l + 5, 1]
    charpoly = np.poly(A)
    assert np.allclose(charpoly, [1, 5, 17, 13])

    def eigenvector(lam):
        return np.array([lam**2 + 5 * lam + 17, lam + 5, 1])

    lam_c, lam_r = -2 + 3j, -1.0
    for lam in (lam_c, lam_r):
        assert np.allclose(A @ eigenvector(lam), lam * eigenvector(lam))
    v = eigenvector(lam_c)
    P = np.column_stack([v.real, v.imag, eigenvector(lam_r)])
    P_inv = np.linalg.inv(P)
    Lam = P_inv @ A @ P
    B_bar, C_bar = P_inv @ B, C @ P
    assert np.allclose(Lam, [[-2, 3, 0], [-3, -2, 0], [0, 0, -1]])
    report(
        "s01a",
        f"characteristic polynomial: {poly_str(charpoly)}\n"
        f"eigenvalues = {show(np.sort_complex(np.linalg.eigvals(A)))}\n"
        f"v(-2 + 3j) = {show(v)}\n"
        f"v(-1) = {show(eigenvector(lam_r))}\n\n"
        f"P =\n{show(P)}\n\n30 P^-1 =\n{show(30 * P_inv)}\n\n"
        f"Lambda = P^-1 A P =\n{show(Lam)}\n\n"
        f"B_bar = P^-1 B =\n{show(B_bar)}\n\nC_bar = C P =\n{show(C_bar)}",
    )

    # b) homogeneous response from a single-mode initial state
    t = np.linspace(0, 3, 601)
    p1, p2, p3 = P.T
    x_osc = np.array([expm(A * tk) @ p1 for tk in t])
    x_exp = np.array([expm(A * tk) @ p3 for tk in t])
    z_osc, z_exp = x_osc @ P_inv.T, x_exp @ P_inv.T
    osc_hand = np.outer(np.exp(-2 * t) * np.cos(3 * t), p1) - np.outer(
        np.exp(-2 * t) * np.sin(3 * t), p2
    )
    exp_hand = np.outer(np.exp(-t), p3)
    report(
        "s01b",
        f"x0 = p1: max |z_3(t)| = {np.abs(z_osc[:, 2]).max():.1e}, "
        f"max |x(t) - e^(-2t)(cos(3t) p1 - sin(3t) p2)| = "
        f"{np.abs(x_osc - osc_hand).max():.1e}\n"
        f"x0 = p3: max |z_1(t)|, |z_2(t)| = {np.abs(z_exp[:, :2]).max():.1e}, "
        f"max |x(t) - e^(-t) p3| = {np.abs(x_exp - exp_hand).max():.1e}",
    )
    assert np.abs(z_osc[:, 2]).max() < 1e-9 and np.abs(z_exp[:, :2]).max() < 1e-9
    assert np.allclose(x_osc, osc_hand) and np.allclose(x_exp, exp_hand)

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 4.6), sharex=True)
    styles = [(BLUE, "-"), (ORANGE, "--"), (AQUA, "-.")]
    cases = [
        (r"$\mathbf{x}_0 = \mathbf{p}_1$ (oscillatory mode)", x_osc, z_osc),
        (r"$\mathbf{x}_0 = \mathbf{p}_3$ (exponential mode)", x_exp, z_exp),
    ]
    for col, (title, x, z) in enumerate(cases):
        for k, (color, ls) in enumerate(styles):
            axes[0, col].plot(t, x[:, k], color=color, ls=ls, label=f"$x_{k + 1}$")
            axes[1, col].plot(t, z[:, k], color=color, ls=ls, label=f"$z_{k + 1}$")
        axes[0, col].set_title(title)
        axes[1, col].set_xlabel("Time")
    axes[0, 0].set_ylabel(r"State $\mathbf{x}(t)$")
    axes[1, 0].set_ylabel(r"Modal coordinates $\mathbf{P}^{-1}\mathbf{x}(t)$")
    for ax in axes.flat:
        ax.grid(True)
        ax.legend(fontsize=8, ncol=3, loc="upper right")
    save(fig, "s01b")

    # c) transfer function matrix from both realizations
    G_orig = channel_tfs(A, B, C)
    G_modal = channel_tfs(Lam, B_bar, C_bar)
    G_full = channel_tfs(A, B, C, cancel=False)
    for key in G_orig:
        assert tf_str(G_orig[key]) == tf_str(G_modal[key])

    def G_hand(s):
        return np.array(
            [[10 / (s + 1), 0], [3 / (s + 1), (s + 4) / (s**2 + 4 * s + 13)]]
        )

    def G_eval(Am, Bm, Cm, s):
        return Cm @ np.linalg.solve(s * np.eye(3) - Am, Bm)

    s_grid = 1j * np.logspace(-2, 2, 200)
    diff_models = max(
        np.abs(G_eval(A, B, C, s) - G_eval(Lam, B_bar, C_bar, s)).max() for s in s_grid
    )
    diff_hand = max(np.abs(G_eval(A, B, C, s) - G_hand(s)).max() for s in s_grid)
    assert diff_models < 1e-9 and diff_hand < 1e-9
    lines = ["original model, C (sI - A)^-1 B before cancellation:"]
    lines += [f"  G{i + 1}{j + 1}(s) = {tf_str(g)}" for (i, j), g in G_full.items()]
    for label, G in (("original", G_orig), ("modal", G_modal)):
        lines.append(f"{label} model after pole-zero cancellation:")
        lines += [f"  G{i + 1}{j + 1}(s) = {tf_str(g)}" for (i, j), g in G.items()]
    lines.append(f"max |G_original(jw) - G_modal(jw)| = {diff_models:.1e}")
    lines.append(f"max |G_original(jw) - G_hand(jw)| = {diff_hand:.1e}")
    report("s01c", "\n".join(lines))

    #########
    # p02 #
    #########
    heading("p02")

    # sample upper block-triangular system: A_1 is 2x2, A_2 is 1x1, A_12 is 2x1
    def A_1(t):
        return np.array([[-1.0, 1 + np.sin(t)], [-0.5, -2 + 0.5 * np.cos(t)]])

    def A_2(t):
        return np.array([[-0.5 + 0.3 * np.cos(t)]])

    def A_12(t):
        return np.array([[1 + 0.5 * t], [np.cos(t)]])

    def A_full(t):
        return np.block([[A_1(t), A_12(t)], [np.zeros((1, 2)), A_2(t)]])

    def fundamental(A_t, n, end):
        sol = solve_ivp(
            lambda s, x: (A_t(s) @ x.reshape(n, n)).ravel(),
            (0, end),
            np.eye(n).ravel(),
            method="DOP853",
            dense_output=True,
            rtol=1e-12,
            atol=1e-13,
        )
        return lambda s: sol.sol(s).reshape(n, n)

    X_1, X_2 = fundamental(A_1, 2, 4.0), fundamental(A_2, 1, 4.0)

    def Phi_1(t, tau):
        return X_1(t) @ np.linalg.inv(X_1(tau))

    def Phi_2(t, tau):
        return X_2(t) @ np.linalg.inv(X_2(tau))

    def Phi_12(t, tau):
        def integrand(s):
            return (Phi_1(t, s) @ A_12(s) @ Phi_2(s, tau)).ravel()

        return quad_vec(integrand, tau, t, epsabs=1e-12, epsrel=1e-11)[0].reshape(2, 1)

    tau = 0.5
    t = np.linspace(tau, 3.5, 31)
    Phi = stm_numeric(A_full, tau, t)
    errors = {
        "Phi_21": max(np.abs(Ph[2:, :2]).max() for Ph in Phi),
        "Phi_11 - Phi_1": max(
            np.abs(Ph[:2, :2] - Phi_1(tk, tau)).max() for tk, Ph in zip(t, Phi)
        ),
        "Phi_22 - Phi_2": max(
            np.abs(Ph[2:, 2:] - Phi_2(tk, tau)).max() for tk, Ph in zip(t, Phi)
        ),
        "Phi_12 - integral formula": max(
            np.abs(Ph[:2, 2:] - Phi_12(tk, tau)).max() for tk, Ph in zip(t, Phi)
        ),
    }
    assert max(errors.values()) < 1e-8
    report(
        "s02",
        "sample system: A_1(t) is 2x2, A_12(t) is 2x1, A_2(t) is 1x1, "
        f"tau = {tau}, t in [{tau}, {t[-1]}]\n"
        + "\n".join(f"max |{k}| = {e:.1e}" for k, e in errors.items()),
    )

    #########
    # p03 #
    #########
    heading("p03")

    def A3(t):
        return np.array([[0.0, t], [0.0, -1.0]])

    def B3(t):
        return np.array([0.0, np.exp(-t)])

    def Phi3(t, tau):
        decay = np.exp(-(t - tau))
        return np.array([[1.0, (tau + 1) - (t + 1) * decay], [0.0, decay]])

    def G3(t, tau):
        return t * ((tau + 1) * np.exp(-tau) - (t + 1) * np.exp(-t)) * (t >= tau)

    def y3(t):
        return t * (3 / np.e - 3 + (2 * np.e * (t + 1) - t**2 - t - 1) * np.exp(-t))

    # a) closed-form STM against the numerical STM
    report(
        "s03a",
        "max relative error of Phi(t, tau) vs numerical STM = "
        f"{stm_error(A3, Phi3):.1e}",
    )

    # b) an impulse at tau sets x(tau+) = B(tau); the output that follows is G(t, tau)
    worst = 0.0
    for tau in (0.0, 1.0, 2.5):
        t = np.linspace(tau, tau + 5, 51)
        sol = solve_ivp(
            lambda s, x: A3(s) @ x,
            (tau, t[-1]),
            B3(tau),
            method="DOP853",
            t_eval=t,
            rtol=1e-11,
            atol=1e-13,
        )
        worst = max(worst, np.abs(t * sol.y[0] - G3(t, tau)).max())
    assert worst < 1e-9
    report(
        "s03b",
        f"max |simulated impulse response - G(t, tau)| for tau = 0, 1, 2.5: {worst:.1e}",
    )

    # c) forced response from x(1) = [1, -2] with u = 1
    x_1 = np.array([1.0, -2.0])
    t = np.linspace(1, 10, 451)
    sol = solve_ivp(
        lambda s, x: A3(s) @ x + B3(s),
        (1, 10),
        x_1,
        method="DOP853",
        t_eval=t,
        rtol=1e-11,
        atol=1e-12,
    )
    y_sim = t * sol.y[0]

    def superposition(t):
        zero_input = t * (Phi3(t, 1.0) @ x_1)[0]
        zero_state = quad(lambda s: G3(t, s), 1, t)[0]
        return zero_input + zero_state

    y_superposition = np.array([superposition(tk) for tk in t])
    assert np.allclose(y3(t), y_sim, atol=1e-8)
    assert np.allclose(y3(t), y_superposition, atol=1e-8)
    report(
        "s03c",
        f"y(1) = {y3(1.0):.6f}\n"
        f"max |y(t) - simulated y(t)| on [1, 10] = {np.abs(y3(t) - y_sim).max():.1e}\n"
        "max |y(t) - (C Phi(t, 1) x(1) + integral of G(t, s) ds)| = "
        f"{np.abs(y3(t) - y_superposition).max():.1e}\n"
        f"y(t) / t -> 3/e - 3 = {3 / np.e - 3:.6f}",
    )

    y_zi = t * (-3 + 2 * (t + 1) * np.exp(1 - t))
    y_zs = t * (3 / np.e - (t**2 + t + 1) * np.exp(-t))
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    ax.plot(t, y3(t), color=BLUE, label=r"$y(t)$, closed form")
    ax.plot(
        t[::15],
        y_sim[::15],
        "o",
        color=BLUE,
        markerfacecolor="white",
        markersize=5,
        label="numerical simulation",
    )
    ax.plot(t, y_zi, color=ORANGE, ls="--", label="zero-input part")
    ax.plot(t, y_zs, color=AQUA, ls="-.", label="zero-state part")
    ax.plot(t, (3 / np.e - 3) * t, color=MUTED, ls=":", label=r"asymptote $(3/e - 3)t$")
    ax.set(xlabel=r"Time $t$", ylabel=r"Output $y(t)$")
    ax.grid(True)
    ax.legend(fontsize=8)
    save(fig, "s03c")

    #########
    # p04 #
    #########
    heading("p04")
    J = np.zeros((5, 5))
    J[:2, :2] = [[-1, 3], [-3, -1]]
    J[2:, 2:] = [[-3, 1, 0], [0, -3, 1], [0, 0, -3]]

    def eJt(t):
        c, s = np.cos(3 * t), np.sin(3 * t)
        E = np.zeros((5, 5))
        E[:2, :2] = np.exp(-t) * np.array([[c, s], [-s, c]])
        E[2:, 2:] = np.exp(-3 * t) * np.array([[1, t, t**2 / 2], [0, 1, t], [0, 0, 1]])
        return E

    # any invertible P works; build an A with the stated structure and check it
    P4 = rng.normal(size=(5, 5))
    P4_inv = np.linalg.inv(P4)
    A4 = P4 @ J @ P4_inv
    r = np.polymul([1, 2, 10], np.poly([-3, -3, -3]))
    ranks = [
        int(
            np.linalg.matrix_rank(
                np.linalg.matrix_power(A4 + 3 * np.eye(5), k), rtol=1e-8
            )
        )
        for k in (1, 2, 3)
    ]
    assert np.allclose(np.poly(A4), r, atol=1e-6) and ranks == [4, 3, 2]
    t = np.linspace(0, 4, 81)
    err_a = max(np.abs(P4_inv @ expm(A4 * tk) @ P4 - eJt(tk)).max() for tk in t)

    def y_b(t):
        return np.array([0, 0, t * np.exp(-3 * t), np.exp(-3 * t), 0])

    def y_c(t):
        return np.array(
            [np.exp(-t) * np.sin(3 * t), np.exp(-t) * np.cos(3 * t), 0, 0, 0]
        )

    err_b = max(np.abs(P4_inv @ expm(A4 * tk) @ P4[:, 3] - y_b(tk)).max() for tk in t)
    err_c = max(np.abs(P4_inv @ expm(A4 * tk) @ P4[:, 1] - y_c(tk)).max() for tk in t)
    assert max(err_a, err_b, err_c) < 1e-9
    report(
        "s04",
        f"r(lambda) = {poly_str(r, 'l')}\n"
        f"charpoly(A) = {poly_str(np.poly(A4), 'l')}  (A = P J P^-1, random P)\n"
        f"rank (A + 3I)^k for k = 1, 2, 3: {ranks}\n"
        f"a) max |P^-1 e^(At) P - closed form| = {err_a:.1e}\n"
        f"b) max |P^-1 e^(At) p_4 - closed form| = {err_b:.1e}\n"
        f"c) max |P^-1 e^(At) p_2 - closed form| = {err_c:.1e}",
    )

    #########
    # p05 #
    #########
    heading("p05")

    # a), b) commuting family: polynomials in one random matrix M
    M = 0.5 * rng.normal(size=(3, 3))
    family = [M, M @ M - M, np.eye(3) + 0.3 * M @ M @ M]
    f = [np.cos, lambda t: t, lambda t: np.exp(-t)]
    F = [np.sin, lambda t: t**2 / 2, lambda t: -np.exp(-t)]  # antiderivatives of f

    def A5(t):
        return sum(f_j(t) * A_j for f_j, A_j in zip(f, family))

    def Phi5_exp(t, tau):
        return expm(sum((F_j(t) - F_j(tau)) * A_j for F_j, A_j in zip(F, family)))

    def Phi5_factored(t, tau):
        factors = [expm((F_j(t) - F_j(tau)) * A_j) for F_j, A_j in zip(F, family)]
        return reduce(np.matmul, factors)

    report(
        "s05b",
        "commuting A_j (polynomials in a random 3x3 M), f = cos(t), t, e^(-t):\n"
        "max relative error vs numerical STM:\n"
        f"  exp(sum_j F_j A_j)       {stm_error(A5, Phi5_exp):.1e}\n"
        f"  prod_j exp(F_j A_j)      {stm_error(A5, Phi5_factored):.1e}",
    )

    # c) A(t) = cos(3t) I + sin(3t) J
    def A5c(t):
        return np.array(
            [[np.cos(3 * t), np.sin(3 * t)], [-np.sin(3 * t), np.cos(3 * t)]]
        )

    def Phi5c(t, tau):
        theta = (np.cos(3 * tau) - np.cos(3 * t)) / 3
        rotation = np.array(
            [[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]]
        )
        return np.exp((np.sin(3 * t) - np.sin(3 * tau)) / 3) * rotation

    samples = rng.uniform(-5, 5, size=(20, 2))
    commutator = max(
        np.abs(A5c(t1) @ A5c(t2) - A5c(t2) @ A5c(t1)).max() for t1, t2 in samples
    )
    report(
        "s05c",
        f"max |A(t1) A(t2) - A(t2) A(t1)| over 20 random pairs = {commutator:.1e}\n"
        "max relative error of Phi(t, tau) vs numerical STM = "
        f"{stm_error(A5c, Phi5c):.1e}",
    )

    #########
    # p06 #
    #########
    heading("p06")

    def A6a(t):
        return np.diag([t, np.cos(t)])

    def Phi6a(t, tau):
        return np.diag([np.exp((t**2 - tau**2) / 2), np.exp(np.sin(t) - np.sin(tau))])

    def A6b(t):
        return np.cos(2 * t) * np.array([[1.0, 1.0], [0.0, 0.0]])

    def Phi6b(t, tau):
        g = np.exp((np.sin(2 * t) - np.sin(2 * tau)) / 2)
        return np.array([[g, g - 1], [0.0, 1.0]])

    def A6c(t):
        return np.array([[-2.0, np.sin(t)], [0.0, -1.0]])

    def Phi6c(t, tau):
        fast, slow = np.exp(-2 * (t - tau)), np.exp(-(t - tau))
        coupling = (
            slow * (np.sin(t) - np.cos(t)) - fast * (np.sin(tau) - np.cos(tau))
        ) / 2
        return np.array([[fast, coupling], [0.0, slow]])

    # exp of the integral of A: not the STM, since A(t) A(s) != A(s) A(t)
    def Phi6c_naive(t, tau):
        return expm(
            np.array([[-2 * (t - tau), np.cos(tau) - np.cos(t)], [0, -(t - tau)]])
        )

    A6d = np.array([[-3.0, -1.0], [-1.0, -3.0]])
    eigvals6, V = np.linalg.eigh(A6d)
    W = np.linalg.inv(V)  # rows are the reciprocal basis

    def Phi6d(t, tau):
        return sum(
            np.exp(lam * (t - tau)) * np.outer(V[:, k], W[k])
            for k, lam in enumerate(eigvals6)
        )

    def Phi6d_hand(t, tau):
        d = t - tau
        return (
            np.exp(-2 * d) * np.array([[1, -1], [-1, 1]])
            + np.exp(-4 * d) * np.array([[1, 1], [1, 1]])
        ) / 2

    commutator = np.abs(A6c(1.0) @ A6c(2.0) - A6c(2.0) @ A6c(1.0)).max()
    naive = max(
        np.abs(Phi6c_naive(t, 0.0) - Phi6c(t, 0.0)).max() for t in np.linspace(0, 3, 31)
    )
    dyadic = max(
        (np.abs(Phi6d(t, 0.0) - expm(A6d * t)) / np.abs(expm(A6d * t)).max()).max()
        for t in np.linspace(-3, 3, 61)
    )
    report(
        "s06",
        "max relative error of Phi(t, tau) vs numerical STM:\n"
        f"  a) {stm_error(A6a, Phi6a):.1e}\n"
        f"  b) {stm_error(A6b, Phi6b):.1e}\n"
        f"  c) {stm_error(A6c, Phi6c):.1e}\n"
        f"  d) {stm_error(lambda t: A6d, Phi6d_hand):.1e}\n"
        f"c) |A(1) A(2) - A(2) A(1)| = {commutator:.3f}\n"
        f"c) max |exp(int A) - Phi| = {naive:.3f}  (exp(int A) is not the STM)\n"
        f"d) eigenvalues = {show(eigvals6)}\n"
        f"d) max relative |dyadic sum - e^(At)| = {dyadic:.1e}",
    )

    #########
    # p07 #
    #########
    heading("p07")

    def A7(t):
        s, c = np.sin(t), np.cos(t)
        return np.array(
            [[-1 + 1.5 * c * c, 1 - 1.5 * s * c], [-1 - 1.5 * s * c, -1 + 1.5 * s * s]]
        )

    def Psi7(t):
        s, c = np.sin(t), np.cos(t)
        return np.array(
            [[np.exp(-t) * s, np.exp(t / 2) * c], [np.exp(-t) * c, -np.exp(t / 2) * s]]
        )

    def dPsi7(t):
        s, c = np.sin(t), np.cos(t)
        return np.array(
            [
                [np.exp(-t) * (c - s), np.exp(t / 2) * (c / 2 - s)],
                [np.exp(-t) * (-s - c), np.exp(t / 2) * (-s / 2 - c)],
            ]
        )

    def Phi7(t, tau):
        return Psi7(t) @ np.linalg.inv(Psi7(tau))

    # a) frozen-time eigenvalues
    t = np.linspace(0, 2 * np.pi, 721)
    frozen = np.array([np.linalg.eigvals(A7(tk)) for tk in t])
    assert np.allclose(frozen.real, -0.25) and np.allclose(
        np.abs(frozen.imag), np.sqrt(7) / 4
    )
    report(
        "s07a",
        f"eigenvalues of A(t) on a 721-point grid over [0, 2 pi]:\n"
        f"  real part in [{frozen.real.min():.6f}, {frozen.real.max():.6f}]\n"
        f"  |imaginary part| in [{np.abs(frozen.imag).min():.6f}, "
        f"{np.abs(frozen.imag).max():.6f}],  sqrt(7)/4 = {np.sqrt(7) / 4:.6f}",
    )

    # b) fundamental solutions
    t = np.linspace(-3, 6, 181)
    residual = max(
        (np.abs(dPsi7(tk) - A7(tk) @ Psi7(tk)) / np.maximum(1, np.abs(Psi7(tk)))).max()
        for tk in t
    )
    det_error = max(abs(np.linalg.det(Psi7(tk)) + np.exp(-tk / 2)) for tk in t)
    assert residual < 1e-12 and det_error < 1e-9
    report(
        "s07b",
        f"max relative |dPsi/dt - A(t) Psi(t)| on [-3, 6] = {residual:.1e}\n"
        f"max |det Psi(t) + e^(-t/2)| on [-3, 6] = {det_error:.1e}",
    )

    # c) STM, monodromy matrix, and divergence
    monodromy = stm_numeric(A7, 0.0, np.array([0.0, np.pi]))[-1]
    multipliers = np.sort(np.linalg.eigvals(monodromy).real)
    assert np.allclose(multipliers, [-np.exp(np.pi / 2), -np.exp(-np.pi)])
    report(
        "s07c",
        "max relative error of Phi(t, tau) vs numerical STM = "
        f"{stm_error(A7, Phi7):.1e}\n"
        f"monodromy matrix Phi(pi, 0) =\n{show(monodromy)}\n"
        f"Floquet multipliers = {show(multipliers)}\n"
        f"-e^(pi/2) = {-np.exp(np.pi / 2):.4f}, -e^(-pi) = {-np.exp(-np.pi):.4f}",
    )

    t = np.linspace(0, 10, 501)
    initial_states = [
        (r"$\mathbf{x}_0 = [1,\ 0]^T = \psi_2(0)$", [1.0, 0.0], BLUE, "-"),
        (r"$\mathbf{x}_0 = [0,\ 1]^T = \psi_1(0)$", [0.0, 1.0], ORANGE, "--"),
        (r"$\mathbf{x}_0 = [1,\ 1]^T/\sqrt{2}$", [2**-0.5, 2**-0.5], AQUA, "-."),
    ]
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    for label, x0, color, ls in initial_states:
        sol = solve_ivp(
            lambda s, x: A7(s) @ x,
            (0, 10),
            x0,
            method="DOP853",
            t_eval=t,
            rtol=1e-11,
            atol=1e-13,
        )
        exact = np.array([Phi7(tk, 0.0) @ x0 for tk in t])
        assert (np.abs(sol.y.T - exact) / np.maximum(1, np.abs(exact))).max() < 1e-7
        ax.semilogy(t, np.linalg.norm(sol.y, axis=0), color=color, ls=ls, label=label)
    ax.semilogy(
        t,
        np.exp(-t / 4),
        color=MUTED,
        ls=":",
        label=r"frozen-eigenvalue decay $e^{-t/4}$",
    )
    ax.set(xlabel=r"Time $t$", ylabel=r"$\|\mathbf{x}(t)\|$")
    ax.grid(True, which="major")
    ax.legend(fontsize=8)
    save(fig, "s07c")


if __name__ == "__main__":
    main()
