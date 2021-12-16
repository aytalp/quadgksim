function Ints = test_func(z)
% Integrands to be tested by the quadgk_sim routine.
% The singularities of the integrands are at the same locations for some of
% the integrands, as in the case of 2-D spectral domain Green's functions.
%
% input:
% z: Complex valued quad points array
% output:
% Ints: Complex valued integrand matrix
% Created by: Aytac Alparslan, 2021.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ints(:,1)=1./((2*z-1)+(z-(0.5-0.1i)));
Ints(:,2)=sin(z)./(2*z-1);
Ints(:,3)=cos(z)./((2*z-1)+(z-0.33));
Ints(:,4)=cos(z).^2./(z-0.33);
Ints(:,5)=cos(z+1)./(2*z-1);
Ints(:,6)=(exp(2i*z)+cos(z))./((2*z-1)+(z-0.33));
Ints(:,7)=(sin(2*z).*cos(z)+3*exp(-1i*z))./((2*z-1).*(z-(0.5-0.1i)));
end