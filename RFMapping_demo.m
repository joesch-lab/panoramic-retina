%% Example script that demonstrates how the receptive fields were calculated
% Note that the script assumes certain parameters that are dependent on the
% experimental setup (eg. stimulus dimensions) and desired resolution
% (eg. latencies)

% As such this script is not expected to run 'out of the box' and and only
% demonstrates an implementation of rf mapping equation in the paper




% this function returns stimulus frame number, tau seconds before time t
sf = precompute_stimulus_pos(stim.times, latencies, response_times); 

%First calculate the 'background' residual RF
mean_response = mean(response, 1); % response is (neurons, time)
filter_bg = zeros(stim.height, stim.width, numel(latencies), 'single'); %(x, y, latency)
for t = 1:size(mean_response,2)
    for tau = 1:numel(latencies)
        stimframe = sf(t, tau);
        if isnan(stimframe); continue; end

        stimulus = single(stim.pattern(:,:,stimframe)); % (x, y, time)
        stimulus = reshape(stimulus, size(stimulus, 1), size(stimulus,2));

        weight = mean_response(t);
        weighted_stim = stimulus*weight;
        filter_bg(:,:,tau) = filter_bg(:,:,tau) + weighted_stim;
    end
end

%Then calculate the 'calcium-triggered-average' for each neuron
filter = zeros(n, stim.height, stim.width, numel(latencies), 'single'); %(neuron, x, y, latency)

for n = 1:size(response,1)
    fprintf("Mapping RF for neuron %d\n", n);
    filter_n = zeros(stim.height, stim.width, numel(latencies), 'single'); 
    
    for t = 1:size(response,2)
        for tau = 1:numel(latencies)
            stimframe = sf(t, tau);
            if isnan(stimframe); continue; end
            
            stimulus = single(stim.pattern(:,:,stimframe)); % (x, y, time)
            stimulus = reshape(stimulus, size(stimulus, 1), size(stimulus,2));
            
            weight = response(n, t);
            weighted_stim = stimulus*weight;
            filter_n(:,:,tau) = filter_n(:,:,tau) + weighted_stim;
        end
    end
    filter(n,:,:,:) = reshape(filter_n, [1, size(filter_n)]); 
end


%Finally, subtract the background
filter_corrected = filter - reshape(filter_bg, [1, size(filter_bg)]); 


function sf = precompute_stimulus_pos(stim.times, latencies, response_times)
%depends on experimental setup and synchronization
%returns nans for the times when there is no data (beginning and end)
return;
end