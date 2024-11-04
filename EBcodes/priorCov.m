classdef priorCov
  % priorCov class
  %  A priorCov object is used to represent a prior covariance matrix Q,
  %       where Q is accessed via function evaluations for both
  %       matrix-vector and matrix-transpose-vector multiplications
  %
  %  The priorCov class has inputs:
  %             xmin, xmax - min and max values for grid
  %             nvec - number of points in the grid
  %             pts - grid points
  %             scale - scalar value
  %             ker_name - name of the kernel (e.g., 'Matern_3by2')
  %             kernel - function handle to evaluate the kernel
  %             gradkernel - 3x1 cell containing function handles to
  %                           evalute the derivatives of the kernel
  %             theta - hyperparameters
  %
  % Authors: Hall-Hooper, Saibaba, Chung, and Miller (2024)
  
    properties
        xmin
        xmax
        nvec
        pts
        scale
        ker_name
        kernel
        gradkernel
        theta
    end
    
    methods
        function Q = priorCov(xmin, xmax, nvec, scale, ker_name)
            Q.xmin = xmin;
            Q.xmax = xmax;
            Q.nvec = nvec;
            Q.scale = scale;
            Q.ker_name = ker_name;
            [Q.kernel,Q.gradkernel] = ker_fcn(ker_name);
            
            [~,Q.pts] = createrow(Q.xmin, Q.xmax, Q.nvec, @(r) exp(-r), Q.scale);
            
            
        end
        
        function [Qm, dQ] = Q_mat(Q, theta)
            Q.theta = theta;
            
            % get dimension of hyperparameters
            r = length(theta);
            
            % stores values of Q matrix for proper use of 'toeplitzproduct' function
            [Qr,~] = createrow(Q.xmin, Q.xmax, Q.nvec, @(r) Q.kernel(r, theta), Q.scale);
            
            % computes the matrix-vector product Qx in an efficient manner; see
            % documentation of toeplitzproduct for details
            Qfun = @(x) toeplitzproduct(x, Qr, Q.nvec);
            
            % stores matrix Q as a funMat object. see documentation of funMat for
            % details
            n = length(Qr);
            Qm = funMat(Qfun, Qfun, [n,n]);
            
            % initialization of cell that holds gradients of Q matrix
            dQ = cell(r,1);
            for j = 1:r
                % stores values of dQ matrix for proper use of 'toeplitzproduct' function
                gradQr = createrow(Q.xmin,Q.xmax,Q.nvec,@(r) Q.gradkernel{j,1}(r, theta),Q.scale);
                
                % computes the matrix-vector product dQx in an efficient manner
                gradQfun = @(x) toeplitzproduct(x, gradQr, Q.nvec);
                
                % stores dQ as a cell of funMat objects
                dQ{j} = funMat(gradQfun, gradQfun, [n,n]);
            end
        end
    
    
        function [U,M] = lowrank(Q, theta, nc)

            sources = Q.pts;
            lims.a = Q.xmin;   lims.b = Q.xmax; 
            compress.method = 'nocompress';
            assert(length(Q.xmin) == 2); % Hard coded for 2D but can be adapted appropriately
            f = @(x1,x2,x3,x4) Q.kernel( sqrt((x3-x1).^2+(x4-x2).^2), theta);
            [U,M] = kernel_lowrank_symm(f, sources, nc, lims, compress);
        end
        
        function [U,M] = lowrank_der(Q, theta, nc)

            sources = Q.pts;
            lims.a = Q.xmin;   lims.b = Q.xmax; 
            compress.method = 'nocompress';
            assert(length(Q.xmin) == 2); % Hard coded for 2D but can be adapted appropriately
            M = cell(1,length(theta));
            for j = 1:length(theta)
                f = @(x1,x2,x3,x4) Q.gradkernel{j}( sqrt((x3-x1).^2+(x4-x2).^2), theta);
                [U,M{j}] = kernel_lowrank_symm(f, sources, nc, lims, compress);
            end
        end
    
    
    end
end