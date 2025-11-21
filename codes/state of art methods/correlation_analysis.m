function correlation_analysis(img, mainTitle)

    % Convert to double
    img = double(img);

    % Extract RGB channels
    R = img(:,:,1);
    G = img(:,:,2);
    B = img(:,:,3);

    figure('Position',[100 100 1400 900]);

    channels = {R, G, B};
    chanNames = {'Red', 'Green', 'Blue'};
    plotIndex = 1;

    % Store correlation values → (3x3 matrix)
    % Rows:   Red, Green, Blue
    % Cols:   Horizontal, Vertical, Diagonal
    results = zeros(3,3);

    for c = 1:3
        channel = channels{c};

        % ---------- HORIZONTAL ----------
        x = channel(:,1:end-1);
        y = channel(:,2:end);
        xh = x(:); yh = y(:);

        subplot(3,3,plotIndex);
        scatter(xh(1:5000), yh(1:5000), 5, chanNames{c}, 'filled');
        title([chanNames{c} ' - Horizontal']);
        xlabel('Pixel(i,j)'), ylabel('Pixel(i,j+1)');
        plotIndex = plotIndex + 1;

        corr_h = corrcoef(xh, yh);
        results(c,1) = corr_h(1,2);


        % ---------- VERTICAL ----------
        x = channel(1:end-1,:);
        y = channel(2:end,:);
        xv = x(:); yv = y(:);

        subplot(3,3,plotIndex);
        scatter(xv(1:5000), yv(1:5000), 5, chanNames{c}, 'filled');
        title([chanNames{c} ' - Vertical']);
        xlabel('Pixel(i,j)'), ylabel('Pixel(i+1,j)');
        plotIndex = plotIndex + 1;

        corr_v = corrcoef(xv, yv);
        results(c,2) = corr_v(1,2);


        % ---------- DIAGONAL ----------
        x = channel(1:end-1,1:end-1);
        y = channel(2:end,2:end);
        xd = x(:); yd = y(:);

        subplot(3,3,plotIndex);
        scatter(xd(1:5000), yd(1:5000), 5, chanNames{c}, 'filled');
        title([chanNames{c} ' - Diagonal']);
        xlabel('Pixel(i,j)'), ylabel('Pixel(i+1,j+1)');
        plotIndex = plotIndex + 1;

        corr_d = corrcoef(xd, yd);
        results(c,3) = corr_d(1,2);
    end

    % Main figure title
    sgtitle(mainTitle, 'FontSize', 16, 'FontWeight', 'bold');

    % ==================== PRINT SUMMARY =====================
    fprintf('\n================= Correlation Summary: %s =================\n', mainTitle);
    fprintf('              Horizontal     Vertical       Diagonal\n');
    fprintf('---------------------------------------------------------------\n');

    fprintf('Red     :     %0.4f        %0.4f         %0.4f\n', ...
            results(1,1), results(1,2), results(1,3));

    fprintf('Green   :     %0.4f        %0.4f         %0.4f\n', ...
            results(2,1), results(2,2), results(2,3));

    fprintf('Blue    :     %0.4f        %0.4f         %0.4f\n', ...
            results(3,1), results(3,2), results(3,3));

    fprintf('===============================================================\n\n');

end
