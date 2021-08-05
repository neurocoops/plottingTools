function [patchX,patchY]=makeErrorBarObject(x,y,e)

uE=y+e;
lE=y-e;

%Make the patch
yP=[lE,fliplr(uE)];
xP=[x,fliplr(x)];

%remove nans otherwise patch won't work
xP(isnan(yP))=[];
yP(isnan(yP))=[];

patchX = xP;
patchY = yP(1,:);
% H.patch=patch(xP,yP(1,:),1,'HandleVisibility','off');