import numpy as np
import os
from numpy.linalg import inv
from scipy.interpolate import interpn
from scipy.optimize import least_squares


def get_unique_numbers(array):
	unique = []
	for number in array:
		if number in unique:
			continue
		else:
			unique.append(number)
	return unique

def non_linear_parameters_95_percent_confidence_interval(fvec, jac):
    """Returns the 95% confidence interval on parameters from
    non-linear fit results."""
    # residual sum of squares
    rss = np.sum(fvec**2)
    # number of data points and parameters
    n, p = jac.shape
    # the statistical degrees of freedom
    nmp = n - p
    # mean residual error
    ssq = rss / nmp
    # the Jacobian
    J = np.matrix(jac)
    # covariance matrix
    c = inv(J.T*J)
    # variance-covariance matrix.
    pcov = c * ssq
    # Diagonal terms provide error estimate based on uncorrelated parameters.
    # The sqrt convert from variance to std. dev. units.
    err = np.sqrt(np.diag(np.abs(pcov))) * 1.96  # std. dev. x 1.96 -> 95% conf
    # Here err is the full 95% area under the normal distribution curve. This
    # means that the plus-minus error is the half of this value
    return err


def norm(cont,obj,theor):
    
    objc_x = []
    objc_y = []
    fluxc_x =[]
    fluxc_y =[]
    objflux = np.zeros((len(obj),2))
     
    for i in range(len(cont)):
        ans = np.where(obj[:,0]==cont[i])
        
        if (np.size(ans)!=0):
            h=20
            s=np.mean(obj[ans[0][0]-h:ans[0][0]+h,1])
            objc_x.append(obj[ans[0][0],0])
            objc_y.append(s)
        
        ans = np.where(theor[:,0]==cont[i])
        
        if (np.size(ans)!=0):
            h=10
            s=max(theor[ans[0][0]-h:ans[0][0]+h,1])
            s1 = np.where(theor[:,1]==s)
            fluxc_x.append(theor[s1[0][0],0])
            fluxc_y.append(s)        
        
    
    objci = np.interp(obj[:,0],objc_x,objc_y)
    fluxci = np.interp(theor[:,0],fluxc_x,fluxc_y)             
    for i in range(len(obj)):
        objflux[i][0] = obj[i][0]
        objflux[i][1] = (obj[i][1]/objci[i])*fluxci[i] 
     
   
    return objflux[:,1]



def funk(x,th,unique_ys,unique_ages,unique_zs,unique_me,obj,err):
    th_interp[:,0]=obj[0:n_elem,0]
    for i in range(n_elem):      
        th_interp[i,1]=(interpn((unique_ys,unique_ages,unique_zs,unique_me),th[:,:,:,:,i],x))
    obj_th_interp=norm(cont,obj[0:n_elem,:],th_interp)
    err_th_interp=norm(cont,err[0:n_elem,:],th_interp)
    print(x)
    f.write('%f %f %f %f\n' % (x[0],x[1],x[2],x[3]))    
                            
    return (obj_th_interp - th_interp[:,1])/err_th_interp 




#----------------- initialisation of parameters 

path_obj = "D:/2023/bta_spec/b514/"  # path to working folder
path_theor = "D:/2022/chimap_prg/spec_mod_arch_5.5A/" # path to a grid of model spectrum

obj=np.loadtxt(path_obj+'Onobjbs_e_r_n.dat') # name of observing spectrum
err=np.loadtxt(path_obj+'Onobj_err_r_n.dat') # name of errors spectrum
cont = np.loadtxt(path_obj+ 'contz001.dat') # list of wavelengths
n_elem = 29000  # number of elements in spectra will use. Or can choose all observing spectra's elements  'len(obj)'



all_theor = os.listdir(path_theor)

ages = []
zs = []
ys = []

#------------ initialisation of grid parameters 
for i in range(len(all_theor)):
    
    age = float(all_theor[i][16:21])
    z = float(all_theor[i][4:10])
    y = float(all_theor[i][11:15])
 
   
    ages.append(age)
    zs.append(z)
    ys.append(y)
    

unique_ages = get_unique_numbers(ages)
unique_zs = get_unique_numbers(zs)
unique_ys = get_unique_numbers(ys)
unique_me = [0,1]

unique_ages.sort()
unique_zs.sort()
unique_ys.sort()



th = np.zeros(shape = (len(unique_ys),len(unique_ages),len(unique_zs),len(unique_me),n_elem))
for i in range(len(th)):
    for j in range(len(th[i])):
        for k in range(len(th[i][j])):
            for q in range(len(th[i][j][k])):
                for l in range(len(th[i][j][k][q])):
                    th[i][j][k][q][l] = -1

ob = np.zeros(shape = (len(unique_ys),len(unique_ages),len(unique_zs),len(unique_me),n_elem))
for i in range(len(ob)):
    for j in range(len(ob[i])):
        for k in range(len(ob[i][j])):
            for q in range(len(ob[i][j][k])):
                for l in range(len(th[i][j][k][q])):
                    th[i][j][k][q][l] = -1

                
er = np.zeros(shape = (len(unique_ys),len(unique_ages),len(unique_zs),len(unique_me),n_elem))
for i in range(len(er)):
    for j in range(len(er[i])):
        for k in range(len(er[i][j])):
            for q in range(len(er[i][j][k])):
                for l in range(len(th[i][j][k][q])):
                    er[i][j][k][q][l] = -1



for u in range(len(all_theor)):
    print('iter=%d' % u)
    theor = np.loadtxt(path_theor + all_theor[u])
    
    age = float(all_theor[u][16:21])
    z = float(all_theor[u][4:10])
    y = float(all_theor[u][11:15])
   
 
    index_age = unique_ages.index(age)
    index_z = unique_zs.index(z)
    index_y = unique_ys.index(y)
    
    if (u%2 == 0):
        index_me=0
    else:
        index_me=1

    
    th[index_y,index_age,index_z,index_me,:] = theor[0:n_elem,1]        
    ob[index_y,index_age,index_z,index_me,:] = norm(cont,obj[:n_elem,:],theor[:n_elem,:])
    er[index_y,index_age,index_z,index_me,:] = norm(cont,err[:n_elem,:],theor[:n_elem,:])
    
print('end')


#----------- finding best fiting from grid of model specrtrum

chi_min=1e15

for i in range(len(ob)):
    for j in range(len(ob[i])):
        for k in range(len(ob[i][j])):
            for q in range(len(ob[i][j][k])):
                
                if (np.sum(th[i,j,k,q,:])==-n_elem):
                    chi=10e15
                else:
                    chi=np.sum(((ob[i,j,k,q,:] - th[i,j,k,q,:])/er[i,j,k,q,:])**2)
            
            
                if (chi<=chi_min):
                    chi_min=chi
                    y_best=unique_ys[i]
                    age_best=unique_ages[j]
                    z_best=unique_zs[k]
                    me_best = unique_me[q]
                    
print('chimin=%f, Y=%f, logAge=%f, Z=%f, me=%f' % (chi_min,y_best,age_best,z_best,me_best))

#----------------- finding of isochron's parameters by non-line minimisation 

th_interp=np.zeros((n_elem,2))
ob_th_interp=np.zeros((n_elem))
err_th_interp=np.zeros((n_elem))

f = open(path_obj+'resultchi_nonline.dat','w') # file with results
f.write('Chimin iso: chimin=%f, Y=%f, logAge=%f, Z=%f, me=%f\n' % (chi_min,y_best,age_best,z_best,me_best))
f.write('---------------------\n')
f.write('Chi nonline\n')
f.write('---------------------\n')


x = np.array([y_best,age_best,z_best,me_best])
lb = np.array([0.23,9.7,0.0001,0])
ub = np.array([0.30,10.15,0.008,1])

res_1 = least_squares(funk,x,loss='linear',bounds=(lb,ub),jac='2-point',verbose=2,args=(th,unique_ys,unique_ages,unique_zs,unique_me,obj,err))


print('Y=%.3f logAge=%.3f Z=%.5f  Me=%.2f' % (res_1.x[0],res_1.x[1],res_1.x[2],res_1.x[3]))

meall=np.array([[-2.5,-2.0,-1.5,-1.0,-1.0,-0.5],[-2.0,-1.5,-1.0,-0.5,-0.5,0]])
fe_h=interpn((unique_me,unique_zs),meall,[res_1.x[3],res_1.x[2]])

print('[Fe/H]=%.2f' % fe_h)

err=non_linear_parameters_95_percent_confidence_interval(res_1.fun, res_1.jac)

print(err)

chi_write = np.sum(res_1.fun**2)

f.write('-----------------------------------------\n')
f.write('Y=%.3f logAge=%.3f Z=%.5f  Me=%.2f chi^2=%.4f\n' % (res_1.x[0],res_1.x[1],res_1.x[2],res_1.x[3],chi_write))
f.write('[Fe/H]=%.2f\n' % fe_h)
f.write('Y_err=%.4f, logAge_err=%.3f, Z_err=%.5f, me_err=%.2f' % (err[0]/2,err[1]/2,err[2]/2,err[3]/2))

f.close()





