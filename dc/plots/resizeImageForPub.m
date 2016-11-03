%  [] = resizeImageForPub(option)
%     option == 'test' || 'onecolumn' || 'potrait' || 'landscape'
function [] = resizeImageForPub(option)

    if strcmpi(option, 'test')
        test_resize;
    end

    hfig = gcf;
    hax = gca;
    funits = hfig.Units;
    aunits = hax.Units;

    hfig.Units = 'inches';
    hfig.Resize = 'off'; % for some reason required.

    hax.Units = 'normalized';
    hax.OuterPosition = [0 0 1 1];

    if strcmpi(option, 'onecolumn')
        hfig.Position(3) = (8.5-1)/2;
    end

    if strcmpi(option, 'portrait')
        hfig.Position(3) = 8.5-1;
    end

    if strcmpi(option, 'landscape')
        hfig.Position(3) = 11.5-1;
    end

    hfig.Units = funits;
    hax.Units = aunits;
    pause(1);hfig.Resize = 'on';
end

function test_resize

    figure; hfig = gcf;
    hfig.Units = 'inches';
    plot(rand([10 10]));
    pbaspect([1.617 1 1]);
    title('asdasd'); xlabel('X'); ylabel('Y');
    beautify;
    drawnow;
    resizeImageForPub('portrait');
    assert(abs(hfig.Position(3) - (8.5-1)) < 10*eps);
    close;
end