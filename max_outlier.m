function output = max_outlier (delta,G_ndiag)
  E = G_ndiag;
  E(G_ndiag>=-delta) = -1.0*G_ndiag(G_ndiag>=-delta);
  E(G_ndiag<-delta) = delta;
  output =  min((sum(G_ndiag+E+speye(size(E)),2)));
end