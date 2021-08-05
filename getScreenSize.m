function [scrSize,width,height]=getScreenSize
scrSize = get(0,'screensize');%use this to get screensize information
% (that way you can adjust the size of your figure depending on your
% screen)
width  = scrSize(3);
height = scrSize(4);