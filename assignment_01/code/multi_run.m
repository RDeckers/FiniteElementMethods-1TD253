 function multi_run(N)
  clf;
  hold on;
  for n = N
    my_first_fem_solver(2^n);
  end;