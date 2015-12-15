function [ eigen ] = problem_02( )
  AK = create_AK_bK([0 0 1], [0 1 0], @f_const, @kappa); %create AK over the reference triangle
  eigen = eig(AK); %and return its eigenvalues
end

