function setwinsize(winhandle,x,y)
% setwinsize(handle,x,y)
%
% Set the size of a window without moving the top left corner
%
% eg: setwinsize(gcf,1000,500)
% Author: from the BICEP/Keck analysis pipeline.

set(winhandle,'Units','pixels');
p=get(winhandle,'Position');
top=p(2)+p(4);
p(3)=x;
p(4)=y;
p(2)=top-p(4);
set(winhandle,'Position',p);
