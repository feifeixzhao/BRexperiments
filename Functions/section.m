function [F,hsrf,hlyr,Z]=section(Z,F,x);

[nz,nx]=size(Z); nt=nz-1;

% by default cell property = time
if(nargin<2 || isempty(F)) F=repmat((1:nt)',[1,nx-1]); end

% by default x = pixels 
if(nargin<3) x=1:nx; end

% plot layers
hsrf=surf(x,Z,0*Z,F); shading flat; view(0,90); grid off;

% plot surfaces
hold on; hlyr=plot(x,Z','k'); hold off; axis tight;

end
