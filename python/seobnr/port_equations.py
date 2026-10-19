"""Regenerate scalar equation ports from the checksum-verified 0.3.7 source.

Run from the repository root, then run Runic on the generated Julia files.
"""
# SPDX-License-Identifier: GPL-3.0-or-later
import ast,re,textwrap,hashlib
from pathlib import Path
import sys
root=Path(sys.argv[1])  # extracted, checksum-verified pySEOBNR 0.3.7 sdist
out=Path('packages/GWSEOBNR/src')
import tomllib
provenance = tomllib.loads(Path("packages/GWSEOBNR/provenance.toml").read_text())
for source, expected in provenance["source_sha256"].items():
    actual = hashlib.sha256((root / source).read_bytes()).hexdigest()
    if actual != expected:
        raise RuntimeError(f"Source checksum mismatch: {source}")
header='# SPDX-License-Identifier: GPL-3.0-or-later\n# Derived from pySEOBNR 0.3.7; see ../NOTICE.md and ../provenance.toml.\n# Modified: generic Julia arithmetic and one-based arrays.\n\n'
class Julia:
 def expr(self,n):
  if isinstance(n,ast.Constant):
   if n.value is True:return 'true'
   if n.value is False:return 'false'
   if n.value is None:return 'nothing'
   if isinstance(n.value,complex):return f'({n.value.real} + {n.value.imag}im)'
   return repr(n.value) if not isinstance(n.value,str) else '"'+n.value+'"'
  if isinstance(n,ast.Name):return {'pi':'π','I':'im','euler_gamma':'EULER_GAMMA','True':'true','False':'false'}.get(n.id,n.id)
  if isinstance(n,ast.Attribute):
   s=ast.unparse(n)
   mp={'self.EOBpars.c_coeffs.a6':'c.a6','self.EOBpars.c_coeffs.dSO':'c.dSO','self.EOBpars.p_params.M':'one(c.nu)','self.EOBpars.p_params.nu':'c.nu','self.EOBpars.p_params.X_1':'c.X_1','self.EOBpars.p_params.X_2':'c.X_2','self._call':'hamiltonian'}
   if s in mp:return mp[s]
   if s.startswith('self.'):return s[5:]
   return s.replace('np.','').replace('fabs','abs')
  if isinstance(n,ast.BinOp):
   op={ast.Add:'+',ast.Sub:'-',ast.Mult:'*',ast.Div:'/',ast.Pow:'^',ast.Mod:'%'}[type(n.op)]
   return f'({self.expr(n.left)} {op} {self.expr(n.right)})'
  if isinstance(n,ast.UnaryOp):return {ast.USub:'-',ast.UAdd:'+',ast.Not:'!'}[type(n.op)]+'('+self.expr(n.operand)+')'
  if isinstance(n,ast.BoolOp):return '('+(' && ' if isinstance(n.op,ast.And) else ' || ').join(map(self.expr,n.values))+')'
  if isinstance(n,ast.Compare):
   op={ast.Eq:'==',ast.NotEq:'!=',ast.Lt:'<',ast.LtE:'<=',ast.Gt:'>',ast.GtE:'>='}
   return '('+self.expr(n.left)+' '+ ' '.join(op[type(o)]+' '+self.expr(v) for o,v in zip(n.ops,n.comparators))+')'
  if isinstance(n,(ast.List,ast.Tuple)):
   a,b=('(',')') if isinstance(n,ast.Tuple) else ('[',']')
   return a+', '.join(map(self.expr,n.elts))+(',' if len(n.elts)==1 and a=='(' else '')+b
  if isinstance(n,ast.Subscript):
   inds=list(n.slice.elts) if isinstance(n.slice,ast.Tuple) else [n.slice]
   val=n.value
   while isinstance(val,ast.Subscript):
    inds=([val.slice] if not isinstance(val.slice,ast.Tuple) else list(val.slice.elts))+inds
    val=val.value
   return self.expr(val)+'['+', '.join(f'{self.expr(i)} + 1' for i in inds)+']'
  if isinstance(n,ast.Call):
   f=self.expr(n.func)
   if f=='evaluate_H':return 'hamiltonian(q, p, c)[1:2]'
   if f=='hamiltonian':return 'hamiltonian(q, p, c)'
   return f+'('+', '.join(map(self.expr,n.args))+')'
  raise ValueError(ast.dump(n))
 def stmts(self,ss,ind=1):
  lines=[];p='    '*ind
  for s in ss:
   if isinstance(s,ast.Expr) and isinstance(s.value,(ast.Constant,ast.Name,ast.Tuple)):continue
   if isinstance(s,ast.Assign):lines.append(p+' = '.join([self.expr(t) for t in s.targets]+[self.expr(s.value)]))
   elif isinstance(s,ast.AugAssign):lines.append(p+self.expr(s.target)+' '+{ast.Add:'+=',ast.Sub:'-=',ast.Mult:'*=',ast.Div:'/='}[type(s.op)]+' '+self.expr(s.value))
   elif isinstance(s,ast.Return):lines.append(p+'return '+self.expr(s.value))
   elif isinstance(s,ast.If):
    cond=self.expr(s.test)
    if isinstance(s.test,(ast.Name,ast.BinOp)):cond='!iszero('+cond+')'
    lines+=[p+'if '+cond]+self.stmts(s.body,ind+1)
    if s.orelse:lines+=[p+'else']+self.stmts(s.orelse,ind+1)
    lines+=[p+'end']
   elif isinstance(s,ast.Raise):lines.append(p+'throw(DomainError(r, "Hamiltonian domain"))')
   elif isinstance(s,ast.Expr):lines.append(p+self.expr(s.value))
   else:raise ValueError(ast.dump(s))
  return lines
j=Julia()
hpath='pyseobnr/eob/hamiltonian/Ham_align_a6_apm_AP15_DP23_gaugeL_Tay_C.pyx'
h=(root/hpath).read_text().splitlines()
functions=[]
for name,a,b in [('hamiltonian',227,339),('hamiltonian_gradient',358,592),('hamiltonian_hessian',610,1284),('hamiltonian_dynamics',1342,1580),('orbital_frequency',1596,1696)]:
 body='\n'.join(h[a:b]); body=body[body.find('        # Coordinate definitions'):] if '# Coordinate definitions' in body else body
 # Start at first declaration after docstring.
 body='\n'.join(h[a:b]); start=next(i for i,l in enumerate(body.splitlines()) if l.strip().startswith('cdef double'))
 body='\n'.join(body.splitlines()[start:]);body=textwrap.dedent(body)
 body=re.sub(r'cdef double\s+','',body)
 if name=='hamiltonian_hessian':
  body=re.sub(r'return np.array\(\[(.*)\]\).reshape\(4, 4\)',r'return (\1)',body)
 tree=ast.parse(body)
 lines=j.stmts(tree.body)
 if name=='hamiltonian_hessian':
  vals=lines[-1].split('return ')[1];lines[-1]='    return SMatrix{4, 4}'+vals
 functions.append('function '+name+'(q, p, c)\n    chi_1, chi_2 = c.chi_1, c.chi_2\n'+ '\n'.join(lines)+'\nend\n')
(out/'hamiltonian.jl').write_text(header+'\n'.join(functions))
fpath='pyseobnr/eob/fits/fits_Hamiltonian.py'
t=ast.parse((root/fpath).read_text());funcs=[]
for f in t.body:
 if not isinstance(f,ast.FunctionDef):continue
 if f.name in ('a6_NS','NR_deltaT_NS'):
  # Convert numpy coefficient array to a regular tuple.
  for s in f.body:
   if isinstance(s,ast.Assign) and isinstance(s.value,ast.Call) and ast.unparse(s.value.func)=='np.array':s.value=s.value.args[0]
 funcs.append('function '+f.name+'('+', '.join(a.arg for a in f.args.args)+')\n'+'\n'.join(j.stmts(f.body))+'\nend\n')
(out/'fits.jl').write_text(header+'\n'.join(funcs)+'''
"""Draw-dependent aligned Hamiltonian coefficients; masses are normalized to total mass."""
function hamiltonian_coefficients(m1, m2, chi1, chi2)
    X_1, X_2, chi_1, chi_2 = promote(m1 / (m1 + m2), m2 / (m1 + m2), chi1, chi2)
    nu = X_1 * X_2
    ap, am = X_1 * chi_1 + X_2 * chi_2, X_1 * chi_1 - X_2 * chi_2
    return (; X_1, X_2, chi_1, chi_2, nu, a6 = a6_NS(nu), dSO = dSO(nu, ap, am))
end
''')
print('generated',[(p.name,len(p.read_text().splitlines())) for p in [out/'fits.jl',out/'hamiltonian.jl']])
# Coefficient routines contain only scalar assignments and piecewise branches.
wpath='pyseobnr/eob/waveform/waveform.pyx'
w=(root/wpath).read_text().splitlines()
coeff=[]
for name,a,b,args in [('rho_coefficients!',410,1145,'nu, dm, a, chiS, chiA, rho_coeffs, rho_coeffs_log, f_coeffs, f_coeffs_vh, extra_PN_terms'),('delta_coefficients!',1170,1310,'nu, dm, a, chiS, chiA, delta_coeffs, delta_coeffs_vh')]:
 # Select from first declaration to next decorator. Strip Cython declarations.
 start=next(i for i in range(a,b) if w[i].strip().startswith('cdef double') and '=' in w[i])
 end=next((i for i in range(start,len(w)) if w[i].startswith('@cython.')),len(w))
 body=textwrap.dedent('\n'.join(w[start:end]))
 body=re.sub(r'cdef double(?: complex)?\s+', '',body)
 tree=ast.parse(body)
 coeff.append('function '+name+'('+args+')\n'+'\n'.join(j.stmts(tree.body))+'\n    return nothing\nend\n')
# Literal angular lookup tables.
allw='\n'.join(w)
for name in ['ylms','ylms_signed']:
 match=re.search(name+r'\[:\] = (\[.*?\n\])',allw,re.S)
 arr=ast.literal_eval(match.group(1))
 coeff.append('const '+name.upper()+' = [\n'+ '\n'.join('    '+' '.join(str(x) for x in row)+';' for row in arr)+'\n]\n')
(out/'flux_coefficients.jl').write_text(header+'const EULER_GAMMA = 0.5772156649015329\n\n'+'\n'.join(coeff))
# Record the source hashes independently of generated Julia formatting.
paths=[hpath,fpath,wpath,'LICENSE','pyseobnr/eob/fits/GSF_fits.py','pyseobnr/eob/fits/IV_fits.py','pyseobnr/eob/fits/MR_fits.py','pyseobnr/eob/utils/nr_utils.py']
Path('packages/GWSEOBNR/provenance.toml').write_text('''upstream = "pyseobnr"
version = "0.3.7"
archive_url = "https://files.pythonhosted.org/packages/dc/ed/7d5b6587f73aa04ffb494bc1a946cbccb9ea659767455606efc89b261525/pyseobnr-0.3.7.tar.gz"
archive_sha256 = "a3125f4b7ba680bb0f0c64f6fc644e1ccb663805d78052b40803f30fd6bbb581"
qnm_version = "0.4.4"

[source_sha256]
'''+''.join('"'+p+'" = "'+hashlib.sha256((root/p).read_bytes()).hexdigest()+'"\n' for p in paths))
# Only the (2,2) output-mode fits are needed by this prototype.
modefits=[]
for filename,prefix in [('IV_fits.py','iv'),('MR_fits.py','mr')]:
 tree=ast.parse((root/'pyseobnr/eob/fits'/filename).read_text())
 cls=next(n for n in tree.body if isinstance(n,ast.ClassDef))
 for f in cls.body:
  if not isinstance(f,ast.FunctionDef) or f.name=='__init__':continue
  expr=next(s.value for s in f.body if isinstance(s,ast.Assign) and isinstance(s.targets[0],ast.Subscript) and ast.unparse(s.targets[0].slice)=='(2, 2)')
  modefits.append('function '+prefix+'_'+f.name+'(c)\n    nu = c.nu\n    dm = c.X_1 - c.X_2\n    chiS = (c.chi_1 + c.chi_2) / 2\n    chiAS = (c.chi_1 - c.chi_2) / 2\n    chi = chiS + chiAS * dm / (1 - 2nu)\n    return '+j.expr(expr)+'\nend\n')
(out/'mode_fits.jl').write_text(header+'\n'.join(modefits))
# Scalar remnant fits. Remove vectorization and input validation (host validates).
tree=ast.parse((root/'pyseobnr/eob/utils/nr_utils.py').read_text())
remnant=[]
for f in tree.body:
 if not isinstance(f,ast.FunctionDef) or f.name not in ['_bbh_HBR2016_ell','bbh_final_spin_non_precessing_HBR2016','bbh_UIBfits_setup','bbh_final_mass_non_precessing_UIB2016']:continue
 ss=[]
 for s in f.body:
  if isinstance(s,ast.Expr) and isinstance(s.value,(ast.Constant,ast.Name,ast.Tuple)):continue
  src=ast.unparse(s)
  if 'np.vectorize' in src or (isinstance(s,ast.If) and src.startswith('if np.any')):continue
  if isinstance(s,ast.Assign) and '_bbh_HBR2016_setup' in src:
   ss+=ast.parse('q = m2 / m1').body;continue
  if isinstance(s,ast.Assign) and src.startswith('sqrt1m4eta ='):
   ss+=ast.parse('sqrt1m4eta = (m1-m2)/(m1+m2)').body;continue
  ss.append(s)
 args=[a.arg for a in f.args.args]
 defaults=[j.expr(d) for d in f.args.defaults]
 for i,d in enumerate(defaults):args[len(args)-len(defaults)+i]+=' = '+d
 remnant.append('function '+f.name+'('+', '.join(args)+')\n'+'\n'.join(j.stmts(ss))+'\nend\n')
(out/'remnant.jl').write_text(header+'''
function calc_isco_radius(a)
    # The exact Schwarzschild limit removes a removable 0/0 in its derivative.
    iszero(primal(a)) && return 6 - 4sqrt(2 / 3) * a
    z1 = 1 + (1 - a^2)^(1 / 3) * ((1 + a)^(1 / 3) + (1 - a)^(1 / 3))
    z2 = sqrt(3a^2 + z1^2)
    return 3 + z2 - sign(a) * sqrt((3 - z1) * (3 + z1 + 2z2))
end

'''+ '\n'.join(remnant))
