import numpy as np
import struct
from lmfit import create_params, fit_report, minimize
from multiprocessing import Pool
import datetime
import sys

start_time = datetime.datetime.now()
print("\nExecution started at: ", start_time)

# def residual(pars, r, data = None, sigma = 1):
#     # for i in range(len(sigma)):
#     #     if sigma[i] == 0:
#     #         sigma[i] = 1
#     vals = pars.valuesdict()
#     A = vals['A']
#     D = vals['D']
#     model = A*r**D
#     if data is None:
#         return model
#     return (model - data)/sigma

def residual(pars, r, L = None, N = 1):
    L = np.where(L == 0.0, 1.0, L)
    vals = pars.valuesdict()
    A = vals['A']
    D = vals['D']
    model = A*r**D
    if L is None:
        return model
    return np.sqrt(N)*(1-model/L)   # tested
    
gs = 0.56
rmin = int(sys.argv[1])
rmax = int(sys.argv[2])
step = int(sys.argv[3])
indir = sys.argv[4]
outfname = f'../fit_results/Nvalues_uniform_hmf_{rmin}-{rmax}.bin'
nthreads = int(sys.argv[5])

r = []
N = []
# L = []
struct_format = 'iiiif'

for i in range(rmin,rmax+1,step):
    r.append(i)
    rval = format(gs*i, '.2f')
    with open(f'{indir}N_values_{rval}Mpc','rb') as f:
        vals = []
        while True:
            data = f.read(struct.calcsize(struct_format))
            if not data:
                break
            integers = struct.unpack(struct_format, data)
            vals.append(integers)
        vals = np.array(vals)
        N.append(vals[:,3])
        # L.append(vals[:,4])

pos = vals[:,:3]
r = np.array(r)
N = np.array(N)
# L = np.array(L)

print('Data read successfully.')

fld = open(outfname, 'wb')

def fit(i):
    # data = L[:,i]
    # for j in range(nr+1):
    #     if N[j,i] == 0:
    #         N[j,i] = 1
    # sigma = L[:,i]/np.sqrt(N[:,i])
    fit_params = create_params(A = 1.0, D = 2.0)
    # out = minimize(residual, fit_params, args=(r,), kws={'data': data, 'sigma': sigma})
    out = minimize(residual, fit_params, args=(r,), 
                    # kws={'L': L[:,i], 'N': N[:,i]}, # for estimating no. of CII emitters
                    kws={'L': N[:,i], 'N': N[:,i]}, # for treating a cell as a unit
                    nan_policy = 'propagate'
                    )
    if out.success != True:
        print(i)
    A = out.params.valuesdict()['A']
    D = out.params.valuesdict()['D']
    reduced = out.redchi
    result = struct.pack('iiifff', int(pos[i,0]), int(pos[i,1]), int(pos[i,2]), A, D, reduced)
    fld.write(result)
    # print(type(residual(out.params, r, L[:,i], N[:,i])))
    # print(f'calculated reduced chi^2 = {np.sum((residual(out.params, r, L[:,i], N[:,i]))**2)/(len(r)-2)}')
    # print(f'x = {pos[i,0]}, y = {pos[i,1]}, z = {pos[i,2]}, A = {A}, D = {D}, reduced chi^2 = {reduced}\n****')

if __name__ == "__main__":
    with Pool(nthreads) as pool:
        pool.map(fit, range(len(N[0])))

fld.close()

end_time = datetime.datetime.now()
print("Execution ended at: ", end_time)
