function done = writeMatrixToFile(fileName, Data)
% Write Lx to "Ld.dat"
h_file = fopen(fileName,'w');
%
[m_L, n_L] = size(Data);
%

fprintf(h_file, '%d\r\n', m_L);

for mm = 1:m_L
    for nn = 1:n_L
%         fprintf(h_file, '%30.16f', Data(mm,nn));
        fprintf(h_file, '%30.16g', Data(mm,nn));
        if mm ~= m_L || nn ~= n_L
            fprintf(h_file, '\t');
        end
    end
    fprintf(h_file, '\r\n');
end
fclose(h_file);
%
done = 1;
end