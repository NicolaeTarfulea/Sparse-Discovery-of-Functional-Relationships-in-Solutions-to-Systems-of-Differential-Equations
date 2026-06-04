function functional_relationships(odesys,t0,X0,T,m,p,S)
% Copyright 2026, All Rights Reserved
% Code by Nicolae Tarfulea
% For paper:"Sparse Discovery of Functional Relationships 
% in Solutions to Systems of Differential Equations"
% by Nicolae Tarfulea
% Input:  odesys is the MATLAB function for the ODE system
%         t0 is the initial value of the independent variable
%         x0 is the initial vector value of the dependent variable at t0
%         T is the end value of the independent variable
%         m is the number of partition intervals of [t0, T]
%         p is the number of Step 3 repetitions
%         S is the number of particular functional relationships 
%           - if they exist- the user wants
% Output: functional relationships and graphical evidence

% Examples of Applications (as presented in Section 4 of the paper) 
% 1. Enzyme Dynamics: 
% functional_relationships(@enzyme,0,[1 0 1 1],1,100,3,6)
% 2. Glycolytic Oscillator: 
% functional_relationships(@glycolytic,0,[1 1 0.1 0.2 0.2 1 0.1],1,10000,3,4)

% The code systematically follows the algorithm outlined in Section 3, 
% executing each step in sequence.

% Step 1: Use the MATLAB built-in function ode45 to compute the numerical 
% solution of the initial value problem associated with the system defined 
% in odesys.
tgrid=linspace(t0,T,m+1); % generates m+1 equally spaced grid points
[t,X]=ode45(odesys,tgrid,X0); % finds the numerical solution

% Enter the candidate functions, namely the augmented library Theta:
LibraryPowers=input(['Do you want to use powers, i.e., X^k? ...' ...
    'Answer 1 for Yes and 0 for No. '])
if LibraryPowers==1
% Enter the nonnegative integer power(s) you want to use as [k1 k2 ...].
    Powers=input('For monomials X^k the powers k = ');
end

% Comment: Likewise, additional functions can be included, such as sin(k*x), 
% cos(k*x), exp(k*x), log(x), etc. For example, the following demonstrates 
% how sin(k*x) can be added for specific values of k:

% LibrarySin=input('Do you want to use sine, i.e., sin(kx)? ...
% Answer 1 for Yes and 0 for No. ')
% if LibrarySin==1
%   SinFreq=input('For sin(kx) the value(s) of k=[k1 k2 ...] are = ')
% end

% Step 2: Generate the numerical matrix Theta_grid.
% First, construct the matrix Theta_grid where each row is obtained 
% by using the numerical solution and the candidate functions. 
[~,c]=size(X); % c is the number of components of the vector solution X
Theta_grid=[]; % creates an empty matrix to store Theta_grid
Tuples=[]; % creates an empty matrix to store Tuples
for k=Powers
  tuples=generateNTuples(k,c); % generates c-tuples with sum k
  [rt,~]=size(tuples);
      for l=1:rt
       C=ones(m+1,1);
         for i=1:c
           C=C.*X(:,i).^tuples(l,i);
         end
       Theta_grid=[Theta_grid, C];
      end
  Tuples=[Tuples;tuples];
end
% Comment: If other functions are considered, e.g., sin(k*x), 
% then Theta_grid must be updated accordingly:
% for k=SinFreq
%    for i=1:c
%       Theta_grid=[Theta_grid sin(k*X(:,i))];
%    end
% end

% Step 3:  
for i=1:p % this for-loop repeats the procedure p times
    % First, each cycle, find the least-squares solutions of Theta_grid*xi=0.
    % That is, solve (Theta_grid'*Theta_grid)*xi=0.
    Trref=rref(Theta_grid'*Theta_grid); % reduced row echellon form
    n=size(Trref,1); % n is the number of rows of Trref
    % The next for-loop detects the rows having just one nonzero (=1) element.
    % The rows with just one nonzero element (=1) correspond to zero
    % xi-coefficients, which are to be discarded.
    zero_xi=[]; % creates an empty vector to store the positions of zero
    % xi-coefficients
    for j=1:n
        k=find(Trref(j,:)~=0);
        if isscalar(k)
            zero_xi=[zero_xi j];
        end
    end
    % Next, discard the xi-coefficients that are zero and upgrade
    % the matrices Theta_grid and Tuples accordingly. This is done
    % by eliminating the columns in Theta_grid and the rows in Tuples
    % that correspond to the xi-coefficients found to be zero.
    Theta_grid(:,zero_xi)=[];
    Tuples(zero_xi,:)=[];
end
% The following if-statement checks whether the matrix Tuples is empty.
% If it is, this corresponds to Step 3, 2(a), in the paper and indicates 
% a failed search.
rows=size(Tuples,1);
if rows==0
    fprintf(['There are no connections of this type. Try a ' ...
        'different library of candidate functions.\n'])
    return
end
% Finally, the code presents the results: the least-squares solutions of
% the system Theta_grid*xi=0, along with the desired relationship F(X)=0.
fprintf(['The xi-coefficients are listed below.\n' ...
    'The xi-coefficients that appear on the right side\n'...
    'of the equations are free parameters;\n' ...
    'the user can choose their values freely.\n'])
GeneralandParticularSolutions(Trref,S,t,Theta_grid) % This function solves 
% and displays the general solution of the linear system Trref*xi=0, 
% togheter with S particular solutions.
fprintf('The relation F(X)=0 is displayed below:\n')
fprintf('xi_1*')
for j=1:c-1
    fprintf('X%d^%d*',j,Tuples(1,j))
end
fprintf('X%d^%d',c,Tuples(1,c))
for i=2:rows
    fprintf('+xi_%d*',i)
    for j=1:c-1
        fprintf('X%d^%d*',j,Tuples(i,j))
    end
    fprintf('X%d^%d',c,Tuples(i,c))
end
% Comment: If the library includes additional functions, such as sin(k*X),
% they should be incorporated into the displayed output accordingly:
% if LibrarySin==1
%     for k=SinFreq
%         for i=1:c
%         rt=rt+1;
%         fprintf('+xi_%d*sin(%d*X%d)',rt,k,i)
%         end
%     end
% end
disp(' = 0')

    function tuples = generateNTuples(n, m)
        % This function generates all possible lists of 
        % m nonnegative integers whose sum is exactly n.

        % First, create an empty matrix to store the tuples
        tuples = [];

        % Call the recursive function to generate the tuples
        generateTuplesRecursively([], n, m);

        % Nested function to generate tuples recursively
        function generateTuplesRecursively(currentTuple, remainingSum, ...
                remainingPositions)
            if remainingPositions == 0
                if remainingSum == 0
                    % If no positions are left and the remaining sum is 0, 
                    % add the tuple
                    tuples = [tuples; currentTuple];
                end
            else
                % Otherwise, iterate over all possible values for the next 
                % position
                for i = 0:remainingSum
                    generateTuplesRecursively([currentTuple, i], ...
                        remainingSum - i, remainingPositions - 1);
                end
            end
        end
        tuples=fliplr(tuples);
    end

    function GeneralandParticularSolutions(A,S,t,B)
    % This function solves the homogeneous system A*xi = 0,
    % displays the general solution, and provides S particular solutions
    % if there are infinitely many solutions.

    [m, n] = size(A);
    R = rref(A); % Reduced row echelon form
    pivotCols = [];

    % Identify pivot columns
    for j = 1:n
        k = find(R(:, j) ~= 0);
        if isscalar(k)
            pivotCols = [pivotCols j];
        end
    end

    % Display general solution
    fprintf('General solution:\n');
    k = 0;
    freeVars = setdiff(1:n, pivotCols);
    for i = pivotCols
        k = k + 1;
        fprintf('xi_%d = ', i);
        for j = i+1:n
            if R(k, j) ~= 0
                fprintf('+(%.4f)*xi_%d ', -R(k, j), j);
            end
        end
        fprintf('\n');
    end
    for j = freeVars
        fprintf('xi_%d = free\n', j);
    end

    % Generate S particular solutions if there are free variables
    N=floor(S/3)+1;
    if ~isempty(freeVars)
        fprintf('\n %d particular solutions:\n',S);
        for s = 1:S
            xi = zeros(n, 1);
            % Assign random integer values to free variables in the range
            % [-5 5], which can be easily updated. 
            % Other pseudo-random number generators can be used, such as 
            % randn(length(freeVars),1)
            randVals = randi([-5 5],length(freeVars), 1);
            xi(freeVars) = randVals;

            % Compute dependent variables using RREF
            for i = 1:length(pivotCols)
                row = i;
                col = pivotCols(i);
                xi(col) = -R(row, freeVars) * xi(freeVars);
            end

            fprintf('Solution %d:\n', s);
            disp(xi');
            subplot(N,2,s) %new
            Cc=B*xi; %new
            plot(t,Cc,'-') %new
            grid on
            title(['Example ', num2str(s)])
            xlabel('t-axis')
            ylabel('F-axis')
        end
        sgtitle('Constraint examples')
    else
        fprintf('\nThe system has only the trivial solution.\n');
    end
end

end