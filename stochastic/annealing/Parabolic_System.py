class Parabolic_System(object):
    
    def __init__(self, q1 ,q2, n, t, N = 32):
        
        import autograd.numpy as np
        from scipy.linalg import expm
        from numpy import concatenate as c


        # Set constants
        self.q1 = q1
        self.q2 = q2
        self.n = n
        self.t = t
        self.N = N

        self.tau = t[1]-t[0]
        self.T = self.t[-1]        



        ## Build constant matrices
        a = np.ones(N,dtype=float)
        d = 4*np.ones(N-1,dtype=float)
        d = c(([2.],d))
        d = c((d,[2.]))
        # d = [2,4,4,...,4,4,2]

        # Linear spline matrix
        # M_{ij} = <psi_i,psi_j>_{L_2(0,1)}
        self.M = np.array( 1./(6.*N) * (np.diag(a, k=-1) + np.diag(d) + np.diag(a, k=1)), dtype=float)
        
        # Stiffness matrix
        # K_{ij} = <psi'_i,psi'_j>_{L_2(0,1)}
        self.K = np.array(np.multiply(N, (np.diag(-a,k=-1) + np.diag(d/2.) + np.diag(-a,k=1))), dtype=float)
        
        L=R=np.zeros((N+1,N+1))
        L[0,0]=1
        R[-1,-1]=1




        ## Compute dynamical system operators      
        self.AN = np.linalg.solve(- self.M,(self.q1*self.K + L),dtype=float)

        self.BN = np.array(np.linalg.solve(self.M,
                                           c(([self.q2],np.zeros(N)))),dtype=float).reshape((N+1,1))
                
        # Discrete-time evolution
        self.ANhat = expm(self.tau*self.AN)

        # Discrete-time input
        self.BNhat= (self.ANhat - np.eye(self.N + 1)) @ np.linalg.solve(self.AN,self.BN)
        
        # Discrete-time output
        self.CNhat = np.zeros(N+1,dtype=float)
        self.CNhat[-1] = 1.
        
        
        
        
        ## Gradient computation
        
        # Derivatives of system operators
        self.dAN_dq1 = np.linalg.solve(-self.M,self.K)        
        self.dBN_dq2 = np.array(np.linalg.solve(self.M,
                                                c(([1],np.zeros(N)))),dtype=float).reshape((N+1,1))

        mm = np.block([[self.AN,self.dAN_dq1],[np.zeros((N+1,N+1)),self.AN]])
        mm = np.multiply(self.tau,mm)
        AdAExp = expm(mm)
        self.ANhat = AdAExp[:(N+1),:(N+1)];
        self.BNhat = (self.ANhat - np.eye(N+1)) @ np.linalg.solve(self.AN,self.BN)
        self.dANhat_dq1 = AdAExp[:(N+1),(N+2):]
        print(self.ANhat.shape, self.BNhat.shape, self.dANhat_dq1.shape)
        
        
        self.dBNhat_dq1 = self.linalg.solve(-self.AN,
                                       self.dAN_dq1 @ 
                                       (self.linalg.solve(self.AN,
                                                          self.ANhat - np.eye(N+1))) - 
                                       self.dANhat_dq1) @ self.BN
        
        self.dBNhat_dq2 = np.linalg.solve(self.AN, (self.ANhat-np.eye(N+1)))*self.dBN_dq2

        
        # State
        X = np.zeros((N+1,n,m+1))
        
        for i in range(m+1):
            X[:,1,i] = np.multiply(self.BNhat , total_u[i,1])
            for j in range(1,n):
                X[:,j,i] = self.ANhat @ X[:,j-1,i] + self.BNhat @ total_u[i,j-1]
        
        
        # Initialize eta for adjoint method
        # goes from t=tau to t=tau*n. At t=0, everything is 0.
        eta = np.zeros[N+1,n+1,m+1]

        for i in range(m+1):
            eta[:,n,i] = (self.CNhat @ X[:,n,i]-Y[i,n]) @ CNhat.T
        
        # Compute gradient contributions
        dJN = np.zeros[1,2+P+1+1]
        
        # Adjoint method
        for j in range(n-1,0,-1):
            for i in range(m+1):
                eta[:,j,i] = self.ANhat.T @ eta[:,j+1,i] + (self.CNhat @ X[:,j,i]-Y[i,j]) @ self.CNhat.T
        
        # Gradient
        for j in range(n+1):
            for i in range(m+1):
                # Gradient of system parameters
                sc = m*(i==m+1)+(i<=m)
                if j>0:
                    for k in range(2):
                        dJN[k] = dJN[k] + sc*(eta[:,j,i].T @ (dAhat_dq(:,:,k)*X(:,j-1,i) + dBhat_dq(:,:,k)*total_u(i,j-1)));
                    end
                end
                JN = JN + sc*(CNhat*X(:,j,i)-Y(i,j))^2;
            end
            
            for r=0:P
                dJN(r+3) = dJN(r+3) + sc*(eta(:,j,m+1)'*Bhat*SplinesP_linear(r+1,j));
            end
        end

        
        
        
        
        
        
        

        """
        # Spline/change of basis matrices for input and state, respectively
        # FIX THIS
        self.phi = np.eye(p,dtype=float)
        
        
        # Build large system operator
        A11 = self.ANhat
        A12 = m(self.BNhat,
                m(c(([1],np.zeros(p-1))).reshape((1,p)), 
                  self.phi))
        A21 = np.zeros((p,self.k(N)),dtype=float)
        A21T = A21.T
        A22 = np.linalg.solve(self.phi,
                              m(np.diag(np.ones(p-1),1),self.phi))
        
        self.A = np.block([[A11,A12],[A21,A22]])
        
        
        # Large spline matrix
        MN = np.block([[self.M, A21T], [A21, self.phi]])
        
        # Build output operator for large system
        self.C = m(c((self.CNhat,np.zeros(p))).reshape((1,p+self.k(N))),
                   np.linalg.inv(MN))
        """
        
    def kern(self, t):
        from numpy.linalg import matrix_power as mp
        if int(t / self.tau) == 0:
            return 0.
        else:
            return self.CNhat @ mp(self.ANhat,int(t/self.tau)-1) @ self.BNhat

    def L(self,i,f):
        import numpy as np
        ti = np.linspace(self.tau, i*self.tau, i)
        if type(f) is np.ndarray:
            fsamp = np.array(f[:i])
        else:
            fsamp = np.array([f(x) for x in ti])
        kernelsamp = np.array([self.kern(x) for x in ti]).reshape(i)
        return np.sum([fsamp[j]*kernelsamp[i-j-1] for j in range(i)])
    
    def fwd(self,u):
        import numpy as np
        N = self.N
        out = np.zeros(self.n)
        s = np.repeat([0],N+1).reshape((N+1,1))
        for i in range(1,self.n):
            s = np.array(self.ANhat @ s + np.multiply(self.BNhat,u[i]))
            out[i] = self.CNhat @ s
        return out
    
    
    
    
    
    
        

    
    
    
    
    
    
    # Dimension of XN
    def k(self,N):
        return N+1
        