function [] = run_gpu(~)

    % % just to initialize GPU for elder versions of matlab
    a = [1,2];
    a = gpuArray(a);
    clear a;
    a = zeros(2,1);
    parfor it = 1:2
        a(it) = it;
    end

end