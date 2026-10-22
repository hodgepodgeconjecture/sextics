"""Minimal sections X=x0+x1 u+x2 u^2, Y=y0+..+y3 u^3 of the cycle-three Jacobian
rational elliptic surface 256Y^2 = 256X^3 + 32 a8 X + b256 (R=1, T=u), msolve format."""
import sys, flint
c1, c2 = map(int, sys.argv[1:3]); char = int(sys.argv[3]); out = sys.argv[4]
names = ['x0','x1','x2','y0','y1','y2','y3','u']
ctx = flint.fmpz_mpoly_ctx.get(tuple(names), 'degrevlex'); g = ctx.gens()
x0,x1,x2,y0,y1,y2,y3,u = g; one = u**0
def a8(c1,c2,al,r,t): return r**2*t*(c1**2*c2**3*(al*r+t)+c1*c2**2*(r+c1**6*t)+c1**3*c2*(r+c2**6*t))
def b256(c1,c2,al,r,t):
    return (-1728*(al*r+t)**2*(r+c1**6*t)**2*(r+c2**6*t)**2
            -32*r**2*t*(c1**3*c2**3*(al*r+t)**2*(r+c2**6*t)+c2**3*(al*r+t)*(r+c1**6*t)**2+c1**3*(r+c1**6*t)*(r+c2**6*t)**2)
            -48*r**2*t*c1**2*c2**2*(al*r+t)*(r+c1**6*t)*(r+c2**6*t)+r**4*t**2*c1**4*c2**4)
X = x0+x1*u+x2*u**2; Y = y0+y1*u+y2*u**2+y3*u**3
E = 256*Y**2-256*X**3-32*a8(c1,c2,1,one,u)*X-b256(c1,c2,1,one,u)
co = {}
for e, cf in E.to_dict().items(): co.setdefault(e[7], {})[e[:7]+(0,)] = cf
polys = [str(ctx.from_dict(co[k])) for k in sorted(co)]
open(out,'w').write(','.join(names[:7])+'\n'+str(char)+'\n'+',\n'.join(polys)+'\n')
print(len(polys), 'equations')
