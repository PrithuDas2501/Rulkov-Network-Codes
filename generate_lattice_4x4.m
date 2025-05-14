function A = generate_lattice_4x4()
    N = 4;
    total = N^2;
    A = zeros(total);

    % Convert 2D indices to 1D
    idx = @(i,j) mod(i-1,N)*N + mod(j-1,N) + 1;

    % Define 8-neighbor offsets
    neighbors = [ -1, -1; -1, 0; -1, 1;
                   0, -1;         0, 1;
                   1, -1;  1, 0;  1, 1 ];

    % Loop through each grid point
    for row = 1:N
        for col = 1:N
            current = idx(row, col);
            for k = 1:size(neighbors,1)
                n_row = mod(row - 1 + neighbors(k,1), N) + 1;
                n_col = mod(col - 1 + neighbors(k,2), N) + 1;
                neighbor = idx(n_row, n_col);
                A(current, neighbor) = 1;
            end
        end
    end
end
