% options - settings struct for the white noise stimulus: The most important ones are:
% options.stimulus_trial_t - Duration of checkers in seconds
% options.update_each_nth_frame - Factor by which to reduce the framerate of display. Eg. 60hz screen, update every 10th frame = 6hz white noise
% options.color_ch - Which color channels to use, Eg. [2,3] uses green and 'blue' leds of the display
% options.specParams.checker_pxsize - Size of each square in pizels (actual visual angle depends on setup)
% options.specParams.random_move - Flag whether to use shifting or stationary checkers
% options.specParams.independentColors- Flag whether each channel has its own pattern and shift for every frame
% options.specParams.move_size_ratio - Factor by which to upsample while shifting. Eg. move_ratio 10 = 150px sized checkers would shift by  multiples of 15px

% The parameters of a previous run are also saved in an 'options' struct in
% the 'targetlogfolder' which can be passed to this function with project = false
% for regnereating the same stimulus for analysis

% stimdata.pattern contains the entire binary pattern of the checkers

function stimdata = ShiftingWN_Demo(options,project,targetlogfolder)% Demo script for the white noise stimulus

if nargin < 3
    targetlogfolder = [pwd  filesep  'stimulii']; %where to save seed etc for regenerating stimulus
end
mkdir(targetlogfolder);


if nargin < 2
    project = true;  %  flag, whether the stimulus is actually shown or just regenerated for analysis
end
debug = true; % flag, whether to show the stimulus fullscreen or not

%%
if project
    stimdata=[];  %nothing to return if projecting
    
    screenid=max(Screen('Screens')); % use highest screen number for presenting
    
    oldskip = Screen('Preference','SkipSyncTests', 0); % 1: continue presentation disregarding sync problems (warning sign will appear!)
    oldVDLevel = Screen('Preference', 'VisualDebugLevel', 1); % make initial screen black
    
    %make offscreen windows and define basic parameters
    if debug %display stim in a smaller window for debugging
        [win , winRect] = Screen('OpenWindow', screenid, 0,[0 0 1280 800],[],[],0);
    else %fullscreen stim for experiments
        [win , winRect] = Screen('OpenWindow', screenid, 0,[],[],[],0);
    end
    
    texture = zeros(winRect(4),winRect(3),3);
    th =size(texture,1);
    tw =size(texture,2);
    rect = [0 0 tw th];
    
    
    texid(1) = Screen('MakeTexture', win, texture, 0, 0, 1);
    texid(2) = Screen('OpenOffscreenWindow', win, [0 0 0],rect);
    ifi = Screen('GetFlipInterval', win);
    vbl=Screen('Flip', win);
    texRect = Screen('Rect', texid(2));
    texRect_size=[texRect(3)-texRect(1) texRect(4)-texRect(2)];
    
else
    if isfield(options,'texRect')
        texRect = options.texRect;
        texRect_size = [texRect(3) texRect(4)];
    else
        texRect=[0,0,1920,1080];
        texRect_size=[1920,1080];
    end
    ifi=options.ifi;
end

if isfield(options,'randomseed') && ~isempty(options.randomseed)
    rng(options.randomseed) %set seed if available
    rseed=options.randomseed;
else
    rseed=rng; %save random seed
end


image_to_present_mono=zeros(texRect_size(2),texRect_size(1),1, 'uint8'); % black image for concatenation
image_to_present=repmat(image_to_present_mono,1,1,3);% RGB black image


drawtex=true; %switch: true=image is made as a matrix and then drawn into texture / false=image components are drawn into texture by PTB commands
subtrials=1; % set number of subtrials



%%


options.specParams.checkers_rows=ceil(texRect_size(2)/options.specParams.checker_pxsize);%calculate number of rows
options.specParams.checkers_columns=ceil(texRect_size(1)/options.specParams.checker_pxsize);%calculate number of rows
oversized_image_to_present=zeros(options.specParams.checker_pxsize.*([options.specParams.checkers_rows options.specParams.checkers_columns])); %initialize oversized image matrix
sizediff=[floor((size(oversized_image_to_present)-size(image_to_present_mono))/2)+[1 1];...
    floor((size(oversized_image_to_present)-size(image_to_present_mono))/2)+size(image_to_present_mono)+[0 0]]; % calculate cropping variables for final image


if ~project %Recreate stimulus for analysis
    if options.specParams.random_move
        rows = options.specParams.move_size_ratio*options.specParams.checkers_rows;
        columns = options.specParams.move_size_ratio*options.specParams.checkers_columns;
        nflips = ceil(size(options.Frameinfo.frame_count,1)/options.update_each_nth_frame);
        stimdata.pattern=zeros(rows,columns,3,nflips, 'logical'); %preallocate
        stimdata.sizediff=options.specParams.move_size_ratio*options.specParams.checker_pxsize;
    else
        rows = options.specParams.checkers_rows;
        columns = options.specParams.checkers_columns;
        nflips = ceil(size(options.Frameinfo.frame_count,1)/options.update_each_nth_frame);
        stimdata.pattern=zeros(rows,columns,3,nflips, 'logical');
        stimdata.sizediff=options.specParams.checker_pxsize;
    end
    framenumber=1;
    stimdata.oversized_image_to_present=oversized_image_to_present;
    stimdata.sizediff=sizediff;
    
    for R=1:options.trials % repeat for R trials %TODO: do this at stim presentation?
        if options.reset_rng_each_trial
            rng(rseed); %reset seed
        end
        while framenumber <= numel(options.Frameinfo.frame_count)
            if mod(framenumber,options.update_each_nth_frame)==1
                flipnumber = (framenumber-1)/options.update_each_nth_frame + 1;
                for c_i=1:numel(options.color_ch)
                    c = options.color_ch(c_i);
                    if c_i==1 || options.specParams.independentColors
                        if isfield(options.specParams, 'random_move') && options.specParams.random_move
                            pattern=randi([0,1],options.specParams.checkers_rows,options.specParams.checkers_columns); %generate pattern matrix
                            pattern = imresize(pattern,size(pattern)*options.specParams.move_size_ratio,'nearest'); %edit
                            shift = randi([0 ,options.specParams.move_size_ratio-1],1,2);
                            pattern = circshift(pattern,shift);
                            stimdata.pattern(:,:,c,flipnumber)=logical(pattern); %generate pattern matrix
                        else
                            stimdata.pattern(:,:,c,flipnumber)=logical(randi([0,1],options.specParams.checkers_rows,options.specParams.checkers_columns,'uint8')); %generate pattern matrix
                        end
                    else
                        c_prev = options.color_ch(c_i-1);
                        stimdata.pattern(:,:,c,flipnumber)= stimdata.pattern(:,:,c_prev,flipnumber);
                    end
                end
            else
            end
            if framenumber==numel(options.Frameinfo.frame_count)
                break
            end
            framenumber=framenumber+1;
        end
    end
    return
end




N_expected_frames=ceil(options.stimulus_trial_t*options.trials/ifi);
waitframes = 1;
frame_count=0;
trialframecount=0;
if project
    Frameinfo.frametime=nan(N_expected_frames,1); %preallocate
    Frameinfo.missed=nan(N_expected_frames,1);
    Frameinfo.frame_count=nan(N_expected_frames,1);
    
    
    if subtrials==1
        fprintf('presenting stimulus %s %s for %d trials with %d seconds each \n',options.type,options.addnote,options.trials,options.stimulus_trial_t)
    else
        fprintf('presenting stimulus %s %s for %d positions/directions and %d trials with %d seconds each \n',options.type,options.addnote,subtrials,options.trials,options.stimulus_trial_t)
    end
    timestart = datetime;
    tsecbegin=GetSecs;
end





for R=1:options.trials % repeat for R trials
    if options.reset_rng_each_trial
        rng(rseed); %reset seed
    end
    new_trial=true;%flag for new trial
    trialframecount=0;
    currentindex = 1; %for movie stimulus
    for Rs=1:subtrials
        new_subtrial=true;%flag for new subtrial
        subtrialframecount = 0;
        if options.wait_intertrial > 0
            Screen('FillRect', texid(2), options.wait_col,texRect);
            Screen('DrawTexture', win, texid(2));
            Screen('Flip', win);
            WaitSecs(options.wait_intertrial);
        end
        
        vbl = GetSecs;
        t0=GetSecs; % reset trial timer
        while GetSecs-t0 < options.stimulus_trial_t % do for time t
            
            if mod(frame_count,options.update_each_nth_frame)==0 %update every nth frame
                if frame_count~=0 && drawtex
                    Screen('Close', presented_texture);
                end
                trialframecount=trialframecount+1;
                for c_i=1:numel(options.color_ch)
                    c = options.color_ch(c_i);
                    if c_i==1 || options.specParams.independentColors
                        pattern=randi([0,1],options.specParams.checkers_rows,options.specParams.checkers_columns); %generate pattern matrix
                        if options.specParams.random_move
                            pattern = imresize(pattern,size(pattern)*options.specParams.move_size_ratio,'nearest'); %edit
                            shift = randi([0 ,options.specParams.move_size_ratio-1],1,2);
                            pattern = circshift(pattern,shift);
                        end
                    else
                        pattern= pattern;
                    end
                    image_to_present(:,:,c)=checker_pattern2image(pattern,c,options,oversized_image_to_present,sizediff);
                end
                if options.specParams.independentColors
                    [r,c]=find((image_to_present(:,:,2)+image_to_present(:,:,3))==options.maxintensity(2));
                    for i = 1:length(r)
                        image_to_present(r(i), c(i),2)=110; % If only green, use 110
                    end
                end
                if drawtex; presented_texture = Screen('MakeTexture', win, image_to_present);end
            end
            % draw the image
            if drawtex; Screen('DrawTexture', texid(2), presented_texture,[],[],0,0);end %draw image to offscreen texture
            Screen('DrawTexture', win, texid(2),texRect,winRect,0,0);%draw image to screen buffer
            
            % finalize and flip on next vertical retrace
            Screen('DrawingFinished', win);
            [vbl,~,~,missed] = Screen('Flip', win, vbl + (waitframes - 0.5) * ifi, 2); % flip next vertical retrace
            
            frame_count=frame_count+1;
            subtrialframecount = subtrialframecount + 1;
            Frameinfo.frame_count(frame_count)=trialframecount;
            Frameinfo.frametime(frame_count)=vbl;
            Frameinfo.missed(frame_count)=missed;
            
        end
    end
end
tsecend=GetSecs;
timeend = datetime;
fprintf("Stimulus presentation over at %f\n", GetSecs)


if drawtex
    Screen('Close', presented_texture);
end
%% save info

if project
    
    fprintf('Saving all data for %s\n', options.type);
    if ~isempty(targetlogfolder)
        outputfile=matfile([targetlogfolder filesep options.type '_' datestr(timestart,'HHMMSS') '.mat'],'writable',true); %create output file
        %save all values
        options.ifi = ifi;
        options.Frameinfo=Frameinfo;
        options.randomseed=rseed;
        options.texRect = texRect;
        
        outputfile.options=options;
        outputfile.timestart=timestart;
        outputfile.tsecbegin = tsecbegin;
        outputfile.timeend=timeend;
        outputfile.tsecend = tsecend;
    end
end
sca;

%%

    function image_to_present=checker_pattern2image(pattern,c,options,oversized_image_to_present,sizediff)
        pattern = pattern.*options.maxintensity(c);
        oversized_image_to_present=imresize(pattern,size(oversized_image_to_present),'nearest');%resize pattern
        image_to_present=oversized_image_to_present(sizediff(1,1):sizediff(2,1),sizediff(1,2):sizediff(2,2));%crop pattern
    end


end