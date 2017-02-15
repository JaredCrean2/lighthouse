
function lighthouse_simple(nu, gamma, omega, t)

  y1 = nu*tan(omega*t)/(gamma - tan(omega*t))
  y2 = gamma*nu*tan(omega*t)/(gamma - tan(omega*t))

  return y1, y2
end

function lighthouse_intermediate(nu, gamma, omega, t)

  v1 = omega*t
  v2 = tan(v1)
  v3 = gamma - v2  # gamma - tan(omega*t)
  v4 = nu*v2  # nu*( tan(omega*t) )
  v5 = v4/v3  # nu*tan(omega*t)/(gamma - v2)
  v6 = v5*gamma
#  v6 = gamma*v5

  y1 = v5
  y2 = v6

  return y1, y2
end

function lighthouse_forward(nu, gamma, omega, t, nu_dot, gamma_dot, omega_dot, t_dot)

  v1 = omega*t
  v1_dot = zero(v1)  # type stability
  v1_dot += t*omega_dot + omega*t_dot

  v2 = tan(v1)
  v2_dot = zero(v2)
  v2_dot += (1/( cos(v1)^2))*v1_dot  # !!! v1_dot

  v3 = gamma - v2
  v3_dot = zero(v3)
  v3_dot += 1*gamma_dot - 1*v2_dot

  v4 = nu*v2
  v4_dot = zero(v4)
  v4_dot += v2*nu_dot + nu*v2_dot

  v5 = v4/v3
  v5_dot = zero(v5)
  v5_dot += (1/v3)*v4_dot + (-v4/(v3^2))*v3_dot
  # this can be simplified using the definition of v5
  # v5_dot = (1/v3)*v4_dot + (-v5/v3)*v3_dot
  # and factoring out 1/v3
  # v5_dot = (v4_dot - v5*v3_dot)/v3

  v6 = gamma*v5
  v6_dot = zero(v6)
  v6_dot += v5*gamma_dot + gamma*v5_dot

  y1 = v5
  y1_dot = v5_dot

  y2 = v6
  y2_dot = v6_dot

  return y1, y2, y1_dot, y2_dot
end


function lighthouse_reverse1(nu, gamma, omega, t, y1_bar, y2_bar)

  # ---------------------------------------------------------------------------
  # forward sweep
  v1 = omega*t
  v2 = tan(v1)
  v3 = gamma - v2  # gamma - tan(omega*t)
  v4 = nu*v2  # nu*( tan(omega*t) )
  v5 = v4/v3  # nu*tan(omega*t)/(gamma - v2)
  v6 = gamma*v5

  y1 = v5
  y2 = v6

  #----------------------------------------------------------------------------
  # reverse sweep

  # initialize all adjoint parts to zero to make accumulation easier
  v6_bar = zero(v6)
  v5_bar = zero(v5)
  v4_bar = zero(v4)
  v3_bar = zero(v3)
  v2_bar = zero(v2)
  v1_bar = zero(v1)
  t_bar = zero(t)
  omega_bar = zero(omega_bar)
  gamma_bar = zero(gamma)
  nu_bar = zero(nu)

  # v6 expression
  # there is a misprint in the paper, v_{-2} should *not* have a bar
  v5_bar += gamma*y2_bar
  gamma_bar += v5*y2_bar

  # v5 expression
  v4_bar += v5_bar(1/v3)
  v3_bar += v5_bar*(-v4/(v3*v3))
  # this can be simplified using the definition of v5
  # v3_bar += v5_bar*(-v5/v3)

  # v4 expression
  nu_bar += v4_bar*v2
  v2_bar += v4_bar*nu

  # v3 expression
  gamma_bar += v3_bar*1
  v2_bar += v3_bar*-1

  # v2 expression
  v1_bar += v2_bar*(1/(cos(v1)^1))

  




  return y1, y2, nu_bar, gamma_bar, omega_bar, t_bar
end
