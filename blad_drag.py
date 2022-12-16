# 
#

l=0.2 # blade length m
rho=1000 # kg/m³
b=0.03 #m
d=0.03 #m
E=3.0*10**8  # pa blade elastic modulus stress/strain
I=b*d**3   # bkade area second moment


U=     # characteristic (wave) velocity scale assumed as magnitude of u at the bed


def Ca(U):
	return rho*b*U**2*I**3/(E*I)