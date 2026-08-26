

from functools import lru_cache
import mpmath as mp

# 设置数值精度（小数位数），比如 50 位
mp.dps = 40



#透反射率

# 右侧输入透射
def right_input_t(k,Delta, g1, g2, phi0, gamma_a = 0):
    """透射振幅 t(k)，使用 mpmath"""

    # phase_shift = mp.e**(mp.j * phi0)
    denom = (k - Delta + 1 * mp.j * gamma_a / 2
             + 2 * mp.j * mp.pi * (abs(g1)**2 + abs(g2)**2)
             + 4 * mp.j * mp.pi * mp.re(g1 * mp.conj(g2)) * mp.e**(mp.j * phi0))
    numer = (k - Delta + 1 * mp.j * gamma_a / 2 - 4 * mp.pi * g1 * mp.conj(g2) * mp.sin(phi0))
    return numer / denom

# 右侧输入透射
def right_input_r(k, Delta, g1, g2, phi0, gamma_a = 0):
    """反射振幅 r(k)，使用 mpmath"""


    # phase_shift = mp.e**(mp.j * phi0)
    denom = (k - Delta + 1 * mp.j * gamma_a / 2
             + 2 * mp.j * mp.pi * (abs(g1)**2 + abs(g2)**2)
             + 4 * mp.j * mp.pi * mp.re(g1 * mp.conj(g2)) * mp.e**(mp.j * phi0) ) 
    numer = (-2 * mp.j * mp.pi *
             (abs(g1)**2
              + abs(g2)**2 * mp.e**(2 * mp.j * phi0)
              + (mp.conj(g1) * g2 + g1 * mp.conj(g2)) * mp.e**(mp.j * phi0) ))
    return numer / denom


# 左侧输入透射
def left_input_r(k, Delta, g1, g2, phi0, gamma_a = 0):
    """
    左入射透射振幅 r_L(p1)：
    r_L = -2π i [|g1|^2 + |g2|^2 e^{-2 i φ0} + (g1* g2 + g1 g2*) e^{-i φ0}] / denom
    其中 denom 与 t_L 相同
    """


    denom = (k - Delta + 1 * mp.j * gamma_a / 2
             + 2 * mp.j * mp.pi * (abs(g1)**2 + abs(g2)**2)
             + 4 * mp.j * mp.pi * mp.re(g1 * mp.conj(g2)) * mp.e**(mp.j * phi0))

    numer = (-2 * mp.j * mp.pi *
             (abs(g1)**2
              + abs(g2)**2 * mp.e**(-2 * mp.j * phi0)
              + (mp.conj(g1) * g2 + g1 * mp.conj(g2)) * mp.e**(-mp.j * phi0)))
    return numer / denom

def left_input_t(k, Delta, g1, g2, phi0, gamma_a = 0):
    """
    左入射透射振幅 t_L(p1)：
    t_L = (p1 - Δ - 4π g1* g2 sin φ0) /
          (p1 - Δ + 2π i (|g1|^2 + |g2|^2) + 4π i Re(g1 g2*) e^{i φ0})
    """
    denom = (k - Delta + 1 * mp.j * gamma_a / 2
             + 2 * mp.j * mp.pi * (abs(g1)**2 + abs(g2)**2)
             + 4 * mp.j * mp.pi * mp.re(g1 * mp.conj(g2)) * mp.e**(mp.j * phi0))
    numer = (k - Delta + 1 * mp.j * gamma_a / 2 - 4 * mp.pi * mp.conj(g1) * g2 * mp.sin(phi0))
    return numer / denom


# 综合函数，根据 direction 选择透射或反射,以及入射方向
def chi(k, Delta, g1, g2, phi0, reflect_transmission_direction='t',input_driection='right', gamma_a = 0):

    if reflect_transmission_direction == 't' and input_driection=='right':
        return right_input_t(k, Delta, g1, g2, phi0, gamma_a)
    elif reflect_transmission_direction == 'r' and input_driection=='right':
        return right_input_r(k, Delta, g1, g2, phi0, gamma_a)    
    elif reflect_transmission_direction == 't' and input_driection=='left':
        return left_input_t(k, Delta, g1, g2, phi0, gamma_a)
    elif reflect_transmission_direction == 'r' and input_driection=='left':
        return left_input_r(k, Delta, g1, g2, phi0, gamma_a)
    else:
        raise ValueError("Invalid direction specified.")
    """综合函数，根据 direction 选择透射或反射"""



# 波形函数
def f_gauss(omega, mu, sigma,T,Np):
    """
    omega：角频率，即被积参数
    mu:均值，平移为0
    sigma:标准差
    T:波之间的时间间隔
    Np:波的数目
    """ 
#     """Normalized Gaussian spectrum f_{μ}(ω)."""
    # return mp.e**(-(omega - mu)**2 / (2 *sigma**2))  / mp.sqrt((sigma * mp.sqrt(mp.pi))) * mp.nsum(lambda m: mp.e**(mp.j*(omega - mu)*m*T), [0, Np-1]) 
    return mp.e**(-(omega)**2 / (2 *sigma**2))  / mp.sqrt((sigma * mp.sqrt(mp.pi))) * mp.nsum(lambda m: mp.e**(mp.j*omega*m*T), [0, Np-1]) 
#     # return mp.e**(-(omega)**2 / (2 *sigma**2))  / mp.sqrt((sigma * mp.sqrt(mp.pi))) * mp.e**(mp.j*omega*T)
#     # return mp.e**(-(omega)**2 / (2 *sigma**2)) / mp.sqrt((sigma * mp.sqrt(mp.pi)))
#     # return mp.e**(-(omega - mu)**2 / (2 *sigma**2)) / mp.sqrt((sigma * mp.sqrt(mp.pi)))



# 一阶关联的计算from functools import lru_cache

def make_psi2(mu, sigma, 
              Delta, g1, g2, phi0, omega0,
            #   reflect_transmission_direction='t',input_driection='right',
              use_infinite_limits=False, L=100, gamma_a = 0):
    """
    Returns a function psi2(t, tau, use_t_left, use_t_right)
    - mu, sigma: Gaussian center and width for f
    - beta, gamma_tot, omega0: parameters for r(ω), t(ω)
    - use_infinite_limits: if True integrate on (-∞,∞); else integrate on [mu-Lσ, mu+Lσ]
    - L: cutoff in units of sigma for finite-range integration
    """
    # mu   = mp.mpf(mu)
    # sigma = mp.mpf(sigma)
    # beta  = mp.mpf(beta)
    # gamma_tot = mp.mpf(gamma_tot)
    # omega0 = mp.mpf(omega0)

    # One-dimensional integral，定义待积分函数
    def I_of_x(x,reflect_transmission_direction,input_driection,T,Np, gamma_a = 0):
        """I(x) = ∫ f(ω) χ(ω) e^{-i ω x} dω.""" #chi(k, Delta, g1, g2, phi0, reflect_transmission_direction='t',inmput_driection='right')
        def integrand(omega):
            return f_gauss(omega , mu, sigma,T,Np) * chi(omega,  Delta, g1, g2, phi0, reflect_transmission_direction,input_driection='right', gamma_a = gamma_a) * mp.e**(-mp.j*(omega)*x)
            # return f_gauss(omega, mu, sigma) * chi_of_omega(omega, use_t, beta, gamma_tot, Gamma0,omega0) * mp.e**(-mp.j*(omega - omega0)*x) 
        
        if use_infinite_limits:
            return mp.quad(integrand, [-mp.inf, mp.inf],maxdegree=10)
            # return mp.quad(integrand, [omega0 - 100, omega0 + 100])
        else:
            a = - L
            b = + L
            return mp.quad(integrand, [a, b])

    # cache I(x, choice) to speed up grids
    @lru_cache(maxsize=None) #通过缓存重复的一维积分来加速
    def I_cached(x_float,T,Np,reflect_transmission_direction,input_driection, gamma_a = 0):
        x = x_float
        return I_of_x(x,reflect_transmission_direction,input_driection,T,Np, gamma_a = gamma_a)
    
    # 最后用上的积分，最后返回的也是这个积分函数
    def psi2(t1, t2, T1, Np1, T2, Np2, test_direction1,test_direction2, gamma_a = 0):
        """
        Compute ψ^{(2)}_{μμ'}(t, t+τ).
        - use_t_left:  False→χ^μ=r，True→χ^μ=t   (for the "left"/ω1 integral)
        - use_t_right: False→χ^{μ'}=r，True→χ^{μ'}=t (for the "right"/ω2 integral)
        """
        x1 = mp.mpf(t1)           # exponent with ω1 reduces to e^{-i ω1 t}
        x2 = mp.mpf(t2)     # exponent with ω2 reduces to e^{-i ω2 (t+τ)}
        phase1 = mp.e**(-mp.j * omega0 * (x1)) #平移项
        phase2 = mp.e**(-mp.j * omega0 * (x2)) #平移项
        I1 = I_cached(x1,T1,Np1,test_direction1,input_driection = 'right', gamma_a = gamma_a) * phase1
        I2 = I_cached(x2,T2,Np2,test_direction2,input_driection = 'right', gamma_a = gamma_a) * phase2
        # 两个积分的时间分别为t1+t2
        
        return (1/(mp.sqrt(2)*mp.pi)) * I1 * I2

    return psi2



#二阶
# def denom_of_right_input_r(k, Delta, g1, g2, phi0):
#     """反射振幅 r(k)，使用 mpmath"""


#     # phase_shift = mp.e**(mp.j * phi0)
#     denom = (k  - Delta
#              + 2 * mp.j * mp.pi * (abs(g1)**2 + abs(g2)**2)
#              + 4 * mp.j * mp.pi * mp.re(g1 * mp.conj(g2)) * mp.e**(mp.j * phi0) + 1e-32) 
#     return -2 * mp.j * mp.pi / denom


def denom_of_right_input_r(k, Delta, g1, g2, phi0, gamma_a = 0):
    """反射振幅 r(k)，使用 mpmath"""


    # phase_shift = mp.e**(mp.j * phi0)
    denom = (k  - Delta + 1 * mp.j * gamma_a / 2
             + 2 * mp.j * mp.pi * (abs(g1)**2 + abs(g2)**2)
             + 4 * mp.j * mp.pi * mp.re(g1 * mp.conj(g2)) * mp.e**(mp.j * phi0) + 1e-32) 
    return -2 * mp.j * mp.pi / denom


@lru_cache(maxsize=None)
def cached_I(x_key, T, Np,
             use_infinite_limits, wmax,
             mu, sigma, Delta, g1, g2, phi0, gamma_a = 0):
    """
    带缓存的一维积分 I(x):
        I(x) = ∫ dω e^{-i ω x} f_gauss(ω) r(ω)
    x_key：用来做 cache key 的实数（比如 float(x) 或 round(...)）
    其余参数保持和原来一致
    """
    x     = mp.mpf(x_key)
    mu    = mp.mpf(mu)
    sigma = mp.mpf(sigma)

    def integrand(omega):
        return (mp.e**(-mp.j * omega * x)
                * f_gauss(omega, mu, sigma, T, Np)
                * denom_of_right_input_r(omega, Delta, g1, g2, phi0, gamma_a))

    if use_infinite_limits:
        I = mp.quad(integrand, [-mp.inf, mp.inf])
    else:
        I = mp.quad(integrand, [-wmax, wmax], maxdegree=10)

    return I 

def N_of_t(t1, t2,
           mu, sigma, Delta, g1, g2, phi0, omega0,T1,Np1,T2,Np2,
           test_direction1, test_direction2,
           use_infinite_limits=True, wmax=50, gamma_a = 0):
    """
    计算 N(t1, t2)。

    参数
    ----
    t1, t2 : float
        时间（与 f_tilde 中 ω 的单位要相容）
    beta, gamma, Gamma0, omega0 : float
        公式中的 β, γ, Γ0, ω0
    f_tilde : callable
        频域包络 \tilde f(ω)，签名为 f_tilde(omega) -> complex (mp.mpf/mp.mpc)
        **建议在函数内使用 mpmath 的 mp 运算，而非 numpy**
    use_infinite_limits : bool
        True: 用 (-∞, +∞) 积分；False: 用有限区间 [-wmax, wmax] 近似
    wmax : float
        有限区间积分的截止频率（当 use_infinite_limits=False 时有效）

    返回
    ----
    complex (mp.mpc)
    """
    # t1 = mp.mpf(t1)
    # t2 = mp.mpf(t2)
    # beta = mp.mpf(beta)
    # gamma = mp.mpf(gamma)
    # Gamma0 = mp.mpf(Gamma0)
    # omega0 = mp.mpf(omega0)

    dt = mp.fabs(t2 - t1)


    # 得写成一个函数
    # def prefactor (g1, g2, phi0,test_direction1, test_direction2):
    # # = -(beta**2) / (2 * mp.pi) * mp.e**(-(gamma/2 + Gamma0) * dt)
    #     Pre_0 = - 1 / (mp.sqrt(2) * mp.pi) * (abs(g1)**2
    #           + abs(g2)**2 * mp.e**(2 * mp.j * phi0)
    #           + (mp.conj(g1) * g2 + g1 * mp.conj(g2)) * mp.e**(mp.j * phi0) )
    #     numer_of_right_input_r =(g1 + g2 * mp.e**(mp.j * phi0) + 1e-16) * (mp.conj(g1) + mp.conj(g2) * mp.e**(mp.j * phi0) + 1e-16)
    #     if test_direction1 == 't' and  test_direction2 == 't':
    #         return Pre_0 * numer_of_right_input_r**2 * (mp.conj(g1) + mp.conj(g2) * mp.e**(-mp.j * phi0) + 1e-16) ** 2/(g1 + g2 * mp.e**(mp.j * phi0) + 1e-16)/(mp.conj(g1) + mp.conj(g2) * mp.e**(mp.j * phi0) + 1e-16)**3
    #     elif test_direction1 == 'r' and  test_direction2 == 'r':
    #         return Pre_0 * numer_of_right_input_r**2 * (mp.conj(g1) + mp.conj(g2) * mp.e**(mp.j * phi0) + 1e-16)/(g1 + g2 * mp.e**(mp.j * phi0) + 1e-16)/(mp.conj(g1) + mp.conj(g2) * mp.e**(mp.j * phi0) + 1e-16)**2
    #     elif test_direction1 == 't' and  test_direction2 == 'r':
    #         return Pre_0 * numer_of_right_input_r**2 * (mp.conj(g1) + mp.conj(g2) * mp.e**(-mp.j * phi0) + 1e-16)/(g1 + g2 * mp.e**(mp.j * phi0) + 1e-16)/(mp.conj(g1) + mp.conj(g2) * mp.e**(mp.j * phi0) + 1e-16)**2
    #     elif test_direction1 == 'r' and  test_direction2 == 't':
    #         return Pre_0 * numer_of_right_input_r**2 * (mp.conj(g1) + mp.conj(g2) * mp.e**(-mp.j * phi0) + 1e-16)/(g1 + g2 * mp.e**(mp.j * phi0) + 1e-16)/(mp.conj(g1) + mp.conj(g2) * mp.e**(mp.j * phi0) + 1e-16)**2
    #     else:
    #         raise ValueError("Invalid direction specified.")

    def prefactor (g1, g2, phi0,test_direction1, test_direction2):
    # = -(beta**2) / (2 * mp.pi) * mp.e**(-(gamma/2 + Gamma0) * dt)
        # Pre_0 = - 1 / (mp.sqrt(2) * mp.pi) * (abs(g1)**2
        #       + abs(g2)**2 * mp.e**(2 * mp.j * phi0)
        #       + (mp.conj(g1) * g2 + g1 * mp.conj(g2)) * mp.e**(mp.j * phi0) )
        Pre_0 = - 1 / (mp.sqrt(2) * mp.pi) 
        # numer_of_right_input_r =(g1 + g2 * mp.e**(mp.j * phi0) + 1e-16) * (mp.conj(g1) + mp.conj(g2) * mp.e**(mp.j * phi0) + 1e-16)
        if test_direction1 == 't' and  test_direction2 == 't':
            return Pre_0 * (mp.conj(g1) + mp.conj(g2) * mp.e**(-mp.j * phi0))**2 * (g1 + g2 * mp.e**(mp.j * phi0))**2 
        elif test_direction1 == 'r' and  test_direction2 == 'r':
            return Pre_0 * (mp.conj(g1) + mp.conj(g2) * mp.e**(mp.j * phi0))**2 * (g1 + g2 * mp.e**(mp.j * phi0))**2
        elif test_direction1 == 't' and  test_direction2 == 'r':
            return Pre_0 * (mp.conj(g1) + mp.conj(g2) * mp.e**(-mp.j * phi0)) * (g1 + g2 * mp.e**(mp.j * phi0))**2 * (mp.conj(g1) + mp.conj(g2) * mp.e**(mp.j * phi0))
        elif test_direction1 == 'r' and  test_direction2 == 't':
            return Pre_0 * (mp.conj(g1) + mp.conj(g2) * mp.e**(-mp.j * phi0)) * (g1 + g2 * mp.e**(mp.j * phi0))**2 * (mp.conj(g1) + mp.conj(g2) * mp.e**(mp.j * phi0))
        else:
            raise ValueError("Invalid direction specified.")
    # decay_term = mp.e**( mp.j * dt * (-Delta + 2 * mp.pi * mp.j * (abs(g1)**2 + abs(g2)**2 + (mp.conj(g1) * g2 + g1 * mp.conj(g2)) * mp.e**(mp.j * phi0) )))
    # def f_tilde(omega,mu,sigma,T,Np):
    #     # f_tilde的高斯函数
    #     mu = mp.mpf(mu)
    #     sigma = mp.mpf(sigma)
    #     norm = 1.0 / mp.sqrt((sigma * mp.sqrt(mp.pi)))
    #     return norm * mp.e**(- (omega)**2 / (2*sigma**2)) * mp.nsum(lambda m: mp.e**(mp.j*omega*m*T), [0, Np-1])
    #     # return norm * mp.e**(- (omega)**2 / (2*sigma**2)) * mp.e**(mp.j*omega*T)
    #     # return norm * mp.e**(- (omega)**2 / (2*sigma**2)) 
    #     # return norm * mp.e**(- (omega -  mu)**2 / (2*sigma**2))  


    # def kernel(omega):
        # (γ/2) / (γ/2 + Γ0 - i(ω - ω0))
        # return (gamma/2) / (gamma/2 + Gamma0 - mp.j*(omega))
        # return 1
        # return (gamma/2) / (gamma/2 + Gamma0 - mp.j*(omega - omega0))
    decay_term = mp.e**(
        mp.j * dt * (
            -Delta + 1 * mp.j * gamma_a / 2
            + 2 * mp.pi * mp.j * (
                abs(g1)**2 + abs(g2)**2
                + (mp.conj(g1) * g2 + g1 * mp.conj(g2)) * mp.e**(mp.j * phi0)
            )
        )
    )

    # ---------- 关键：只对 x = (t1+t2-dt)/2 做一次积分，并缓存 ----------
    # (t1 + t2 - dt)/2 = min(t1, t2) 这一点之前我们已经推过
    x = (t1 + t2 - dt) / 2
    # 为了避免浮点误差导致 key 不一致，可以适当 round 一下
    x_key = float(round(x, 12))

    # 整体相位 e^{-i ω0 (t1+t2)/2}（两套代码已经统一）
    phase1 = mp.e**(-mp.j * omega0 * (t1 + t2) / 2)

    # 调用带缓存的一维积分
    I1 = cached_I(x_key, T1, Np1,
                  use_infinite_limits, wmax,
                  mu, sigma, Delta, g1, g2, phi0, gamma_a) * phase1
    return prefactor(g1, g2, phi0,test_direction1, test_direction2) * decay_term * (I1 * I1) 



# 测试代码
if __name__ == "__main__":
    k = 1
    Delta = 0.5
    g1 = 0.1
    g2 = 0.1
    phi0 = mp.pi / 4
    gamma_a = 0.01
    reflect_transmission_direction = 'r'  # 't' for transmission, 'r' for reflection
    input_driection = 'right'  # 'right' for right input, 'left' for left input
    print(chi(k, Delta, g1, g2, phi0, reflect_transmission_direction, input_driection, gamma_a))
    print(chi(k, Delta, g1, g2, phi0, reflect_transmission_direction, input_driection, gamma_a = 0))

    k = 5.8
    Delta = -0.01
    g1 = 0.0282
    g2 = 0.0282
    phi0 = mp.pi / 2
    gamma_a = 0.0
    test_direction1 = 'r'  # 't' for transmission, 'r' for reflection
    test_direction2 = 'r'  # 't' for transmission, 'r' for reflection
    input_driection = 'right'  # 'right' for right input, 'left' for left input
    sigma = 0.01
    omega0= 5.8
    t1 = 0.8
    t2 = 0.8
    delta_t = 1e2
    T1 = 1
    T2 = 1
    Np1 = 0
    Np2 = 0
    #     psi2 = make_psi2(mu, sigma,
    #                      Delta, g1, g2, phi0, omega0,
    #                      use_infinite_limits=use_infinite_limits, L=wmax)
    psi2_test1 = make_psi2(k, sigma,
                    Delta, g1, g2, phi0, omega0,
                    use_infinite_limits = "true", L=40, gamma_a=gamma_a)
    Nval1_test1 = psi2_test1(t1 * delta_t, t2 * delta_t, T1, Np1, T2, Np2,
                            test_direction1, test_direction2, gamma_a=gamma_a)

    psi2_test2 = make_psi2(k, sigma,
                    Delta, g1, g2, phi0, omega0,
                    use_infinite_limits = "true", L=40, gamma_a=0)
    Nval1_test2 = psi2_test2(t1 * delta_t, t2 * delta_t, T1, Np1, T2, Np2,
                            test_direction1, test_direction2,gamma_a=0)


    print(abs(Nval1_test1) ** 2 * 1e6,abs(Nval1_test2) ** 2 * 1e6)



    Nval2_test1 = N_of_t(t1 * delta_t, t2 * delta_t, k, sigma,
                           Delta, g1, g2, phi0, omega0, T1, Np1, T2, Np2,
                           test_direction1, test_direction2,
                           use_infinite_limits="true", wmax=40, gamma_a=gamma_a)
    Nval2_test2 = N_of_t(t1 * delta_t, t2 * delta_t, k, sigma,
                            Delta, g1, g2, phi0, omega0, T1, Np1, T2, Np2,
                            test_direction1, test_direction2,
                            use_infinite_limits="true", wmax=40, gamma_a=0)

    print(abs(Nval2_test1) ** 2 * 1e6,abs(Nval2_test2) ** 2 * 1e6)
    print(abs(Nval1_test1 + Nval2_test1) ** 2 * 1e6 , abs(Nval2_test1 + Nval2_test2) ** 2 * 1e6)



















