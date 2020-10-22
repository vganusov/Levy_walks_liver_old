function h=tailplot(x, xmin, alpha)
% 

% reshape input vector
x = reshape(x,numel(x),1);
% initialize storage for output handles
h = zeros(2,1);

% select method (discrete or continuous) for plotting
if     isempty(setdiff(x,floor(x))), f_dattype = 'INTS';
elseif isreal(x),    f_dattype = 'REAL';
else                 f_dattype = 'UNKN'; 
end;
if strcmp(f_dattype,'INTS') && min(x) > 50,
    f_dattype = 'REAL';
end;

% estimate xmin and alpha, accordingly
switch f_dattype,
    
    case 'REAL',
        n = length(x);
        c = [sort(x) (n:-1:1)'./n];
        q = sort(x(x>=xmin));
        cf = [q (q./xmin).^(1-alpha)];
        cf(:,2) = cf(:,2) .* c(find(c(:,1)>=xmin,1,'first'),2);

        figure(1);subplot(1,4,2)
        h(1) = loglog(c(:,1),c(:,2),'bo','MarkerSize',8,'MarkerFaceColor',[1 1 1]); hold on;
        h(2) = loglog(cf(:,1),cf(:,2),'k--','LineWidth',2); hold off;
        xr  = [10.^floor(log10(min(x))) 10.^ceil(log10(max(x)))];
        xrt = (round(log10(xr(1))):2:round(log10(xr(2))));
        if length(xrt)<4, xrt = (round(log10(xr(1))):1:round(log10(xr(2)))); end;
        yr  = [10.^floor(log10(1/n)) 1];
        yrt = (round(log10(yr(1))):2:round(log10(yr(2))));
        if length(yrt)<4, yrt = (round(log10(yr(1))):1:round(log10(yr(2)))); end;
        set(gca,'XLim',xr,'XTick',10.^xrt);
        set(gca,'YLim',yr,'YTick',10.^yrt);%,'FontSize',9.35);
        ylabel('Prb(D \geq d)');%,'FontSize',9.35);
        xlabel('Branch lengths d (\mum)');%,'FontSize',9.35)

    case 'INTS', 
        n = length(x);        
        q = unique(x);
        c = hist(x,q)'./n;
        c = [[q; q(end)+1] 1-[0; cumsum(c)]]; c(c(:,2)<10^-10,:) = [];
        cf = ((xmin:q(end))'.^-alpha)./(zeta(alpha) - sum((1:xmin-1).^-alpha));
        cf = [(xmin:q(end)+1)' 1-[0; cumsum(cf)]];
        cf(:,2) = cf(:,2) .* c(c(:,1)==xmin,2);

        figure();
        h(1) = loglog(c(:,1),c(:,2),'bo','MarkerSize',8,'MarkerFaceColor',[1 1 1]); hold on;
        h(2) = loglog(cf(:,1),cf(:,2),'k--','LineWidth',2); hold off;
        xr  = [10.^floor(log10(min(x))) 10.^ceil(log10(max(x)))];
        xrt = (round(log10(xr(1))):2:round(log10(xr(2))));
        if length(xrt)<4, xrt = (round(log10(xr(1))):1:round(log10(xr(2)))); end;
        yr  = [10.^floor(log10(1/n)) 1];
        yrt = (round(log10(yr(1))):2:round(log10(yr(2))));
        if length(yrt)<4, yrt = (round(log10(yr(1))):1:round(log10(yr(2)))); end;
        set(gca,'XLim',xr,'XTick',10.^xrt);
        set(gca,'YLim',yr,'YTick',10.^yrt,'FontSize',16);
        ylabel('Pr(X \geq x)','FontSize',16);
        xlabel('x','FontSize',16)

    otherwise,
        fprintf('(PLPLOT) Error: x must contain only reals or only integers.\n');
        h = [];
        return;
end;
