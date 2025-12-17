function plot_hist(img, titleStr)

    figure('Position',[200 200 1200 400]);   % wide figure for 3 subplots

    chanColors = {'r','g','b'};              % plot colors
    chanNames  = {'Red','Green','Blue'};     % titles

    for i = 1:3
        subplot(1,3,i);

        [counts, bins] = imhist(img(:,:,i));

        h = stem(bins, counts, 'LineWidth', 1.1, 'Color', chanColors{i});
        set(h, 'Marker', 'none');       % <-- removes round circles

        title([titleStr ' - ' chanNames{i}]);
        xlabel('Pixel value');
        ylabel('Count');

        axis tight;                          
        ylim([0 max(counts)*1.1]);      % avoid cutoff
        grid on;
    end
end
