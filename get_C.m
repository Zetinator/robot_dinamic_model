function C = get_C(D,q,q_d)
  n_links = size(D,1);

  aux_sum = 0;
  for k=1:n_links
    for j=1:n_links
      for i=1:n_links
        c(i,j,k) = (1/2)*(diff(D(k,j),q(i))+diff(D(k,i),q(j))-diff(D(i,j),q(k)))*q_d(i);
      end
      C(k,j) = sum(c(:,j,k));
    end
  end

end
