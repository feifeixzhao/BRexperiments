function [S,Z]=plot_Section(Z,x,smpInterval);

%topography 2 stratigraphy
S = topo2strat(Z);

%plot a section, colormaped by depositional time
[F,hsrf,hlyr,Z]=section(S(1:smpInterval:end,:),[],x);