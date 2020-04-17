function draw_bistatic_ellipse_overlap(xtx,xrx,xtgt,dist_err)
% Draw a set of bistatic isorange ellipses, given a transmitter at xtx, set
% of receivers at xrx, and target at xtgt, with thickness provided by
% dist_err.
%
% A separate set of ellipses will be drawn for each target.
%
% Nicholas O'Donoughue
% 24 Feb 2020

figure;
plot(xtx(1),xtx(2),'+','DisplayName','Transmitter');
hold on;
plot(xrx(1,:),xrx(2,:),'+','DisplayName','Recievers');
plot(xtgt(1,:),xtgt(2,:),'^','DisplayName','Targets');

for idx_tgt = 1:size(xtgt,2)
    this_tgt = xtgt(:,idx_tgt);
    
    for idx_rx = 1:size(xrx,2)
        this_rx = xrx(:,idx_rx);
        
        % Plot true position ellipse
        ell = utils.drawEllipseFromFoci([xtx,this_rx],this_tgt,[],0);
        hdl=plot(ell(1,:),ell(2,:),'k-','DisplayName','True Solution');
        if idx_tgt ~= 1 || idx_rx ~=1
            hdl.HandleVisibility = 'off';
        end
        
        % Plot error ellipse
        ell_pos = utils.drawEllipseFromFoci([xtx,this_rx],this_tgt,[],dist_err);
        ell_neg = utils.drawEllipseFromFoci([xtx,this_rx],this_tgt,[],-dist_err);
        hdl=plot(ell_pos(1,:),ell_pos(2,:),'k-.','DisplayName','Error Bounds');
        if idx_tgt ~=1 || idx_rx ~=1
            hdl.HandleVisibility = 'off';
        end
        hdl=plot(ell_neg(1,:),ell_neg(2,:),'k-.','DisplayName','Error Bounds');
        hdl.HandleVisibility = 'off';
    end
end