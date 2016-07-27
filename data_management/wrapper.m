function output_matrix = wrapper(input_matrix)
	last_nonzero_colum = max(find(max(input_matrix)));
	output_matrix = input_matrix(:,1:last_nonzero_colum);
end
