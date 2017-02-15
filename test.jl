using FactCheck

include("lighthouse.jl")

# sample values to test
nu_vals = linspace(0.25, 1.25, 4)
gamma_vals = linspace(0.5, 3, 4)
omega_vals = linspace(1, 3, 4)
t_vals = linspace(0, 5, 4)

facts("----- Testing lighthouse function -----") do

for nu in nu_vals
  for gamma in gamma_vals
    for omega in omega_vals
      for t in t_vals
        y1_simple, y2_simple = lighthouse_simple(nu, gamma, omega, t)
        y1_intermediate, y2_intermediate = lighthouse_intermediate(nu, gamma, omega, t)
        @fact y1_simple --> roughly(y1_intermediate, atol=1e-13)
        @fact y2_simple --> roughly(y2_intermediate, atol=1e-13)
      end
    end
  end
end

end


facts("----- Testing lighthouse forward mode -----") do
h = 1e-20
pert = Complex128(0, h)
for nu in nu_vals
  for gamma in gamma_vals
    for omega in omega_vals
      for t in t_vals

        # test first variable
        y1_cs, y2_cs = lighthouse_intermediate(nu + pert, gamma, omega, t)
        y1_fm, y2_fm, y1_fmdot, y2_fmdot = lighthouse_forward(nu, gamma, omega, t, 1, 0, 0, 0)

        @fact real(y1_cs) --> roughly(y1_fm, atol=1e-13)
        @fact real(y2_cs) --> roughly(y2_fm, atol=1e-13)
        @fact imag(y1_cs)/h --> roughly(y1_fmdot, atol=1e-13)
        @fact imag(y2_cs)/h --> roughly(y2_fmdot, atol=1e-13)

        # test second variable
        y1_cs, y2_cs = lighthouse_intermediate(nu, gamma + pert, omega, t)
        y1_fm, y2_fm, y1_fmdot, y2_fmdot = lighthouse_forward(nu, gamma, omega, t, 0, 1, 0, 0)

        @fact real(y1_cs) --> roughly(y1_fm, atol=1e-13)
        @fact real(y2_cs) --> roughly(y2_fm, atol=1e-13)
        @fact imag(y1_cs)/h --> roughly(y1_fmdot, atol=1e-13)
        @fact imag(y2_cs)/h --> roughly(y2_fmdot, atol=1e-13)

        # test 3rd variable
        y1_cs, y2_cs = lighthouse_intermediate(nu, gamma, omega + pert, t)
        y1_fm, y2_fm, y1_fmdot, y2_fmdot = lighthouse_forward(nu, gamma, omega, t, 0, 0, 1, 0)

        @fact real(y1_cs) --> roughly(y1_fm, atol=1e-13)
        @fact real(y2_cs) --> roughly(y2_fm, atol=1e-13)
        @fact imag(y1_cs)/h --> roughly(y1_fmdot, atol=1e-13)
        @fact imag(y2_cs)/h --> roughly(y2_fmdot, atol=1e-13)

        # test 4th variable
        y1_cs, y2_cs = lighthouse_intermediate(nu, gamma, omega, t + pert)
        y1_fm, y2_fm, y1_fmdot, y2_fmdot = lighthouse_forward(nu, gamma, omega, t, 0, 0, 0, 1)

        @fact real(y1_cs) --> roughly(y1_fm, atol=1e-13)
        @fact real(y2_cs) --> roughly(y2_fm, atol=1e-13)
        @fact imag(y1_cs)/h --> roughly(y1_fmdot, atol=1e-13)
        @fact imag(y2_cs)/h --> roughly(y2_fmdot, atol=1e-13)

      end
    end
  end
end

end


facts("----- Testing lighthouse reverse mode -----") do
h = 1e-20
pert = Complex128(0, h)

jac_forward = zeros(2, 4)
jac_reverse = zeros(2, 4)
for nu in nu_vals
  for gamma in gamma_vals
    for omega in omega_vals
      for t in t_vals
        # compute the jacobian

        y1, y2, nu_bar, gamma_bar, omega_bar, t_bar = lighthouse_reverse1(nu, gamma, omega, t, 1, 0)





      end
    end
  end
end

end
# construct jacobian by rows, columns
