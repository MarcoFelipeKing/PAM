function trans = countFcn(seqs)
    %# accumulate transition matrix from all sequences
    trans = zeros(5,5);
    for i=1:numel(seqs)
        trans = trans + sparse(seqs{i}(1:end-1), seqs{i}(2:end), 1, 5,5);
    end

    %# normalize into proper probabilities
    trans = bsxfun(@rdivide, trans, sum(trans,2));
end
