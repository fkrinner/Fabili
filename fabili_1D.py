from math import sin, cos

class fabili1D:
	def __init__(self, funct):
		self.funct    = funct
		self.stepSize = .001
		self.numLim   = 1.e-6

	def getEstimate(self, x):
		f   = self.funct.f(x)
		Df  = self.funct.Df(x)
		DDf = self.funct.DDf(x)
		A   = DDf/2
		B   = Df - DDf * x
		if self.isZero(DDf):
			if self.isZero(Df):
				raise ValueError("No derivative in any direction")
			return x + self.stepSize * Df/abs(Df)
		if DDf > 0.:
			return -B/(2*A)
		else:
			BB = 4*A*x + B
			return BB/(2*A) # AA = -A		

	def isZero(self, v):
		if abs(v) < self.numLim:
			return True
		return False


class function:
	def __init__(self, a,b,c):
		self.A = a
		self.B = b
		self.C = c

	def f(self, x):
		return self.A*x**2 + self.B*x + self.C
	
	def Df(self, x):
		return 2*self.A*x + self.B

	def DDf(self, x):
		return 2*self.A

class sine:
	def __init__(self):
		pass

	def f(self, x):
		return sin(x)

	def Df(self, x):
		return cos(x)

	def DDf(self, x):
		return -sin(x)

class cose:
	def __init__(self):
		pass

	def f(self, x):
		return cos(x)

	def Df(self, x):
		return -sin(x)

	def DDf(self, x):
		return -cos(x)


def main():
	ff = cose()	
	x0 = 0.001
	fab = fabili1D(ff)
	nn = 20
	for _ in range(nn):
		x  = fab.getEstimate(x0)
		print x,ff.Df(x)
		x0 = x
	print ff.f(x0), ff.Df(x0)

if __name__ == "__main__":
	main()
	
