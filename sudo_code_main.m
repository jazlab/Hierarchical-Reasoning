function [] = sudo_code_main();
%% THe purpose of this sudo code.
% in this code, we will show how to use different sub-functions to fit the
% models. It is written as structured as possible to show how to use
% different sub-function to fit models to the behavior.

addpath('/Users/mortezaimac/Dropbox (MIT)/PhD-MIT-BCS/Experiments/Codes/Electrophysiology/Paper1/');
close all; clc;


test_plot_enable = 1; % enable test plots
%% datapath, subject info, ...

Subject.Names = {'Kodaly', 'Ives'};
FullInternal_Enable = 0; % 0: rule report 1: no rule report
stim_switch_model_enable = 0;

% All Kodaly's behavioral data (during ephys and learned data):
Subject.SessionDate(1).Date = { '20170306', '20170310-2', '20170310', '20170312', '20170314', '20170315', '20170316', '20170317', '20170320', '20170321', '20170322', '20170630', '20170701', '20170703', '20170714', '20170715', '20170716', '20170717', '20170721', '20170723', '20170724', '20170726', '20170727' , '20170731', '20170801', '20170802', '20170803', '20170804', '20170806', '20170808', '20170809', '20170817', '20170819', '20170821', '20170825', '20170828', '20170829', '20170831', '20170901', '20170911', '20170913', '20170914', '20170915' };
% All Ives's behavioral data (during ephys and learned data):
Subject.SessionDate(2).Date = {'20170517', '20170831', '20170921', '20170328', '20170524', '20171101', '20171107', '20180104', '20180109', '20180205', '20170829', '20170522', '20170520', '20170520', '20170307-2', '20170410', '20170922', '20171012', '20171021', '20171023',  '20171117', '20171119', '20180101', '20180308', '20180314', '20180327', '20180329'};


myDir = '/Users/mortezaimac/Dropbox (MIT)/PhD-MIT-BCS/Experiments/';
myProcessedDataFolder = 'MonkeyPsychData/Behavior_ExtractedData/';

iSubject = 1; % subject's ID
iDate = 1; % session number
disp(strcat('[iSubject: ', num2str(iSubject), '---iDate:  ', num2str(iDate), ']'));



%% loading dataset (Here it is one session)
myExtractedDataName = strcat(myDir, myProcessedDataFolder, Subject.Names{iSubject}, '-', Subject.SessionDate(iSubject).Date{iDate},'.mat');
load(myExtractedDataName);

InputAll(iSubject, iDate).StateExp = double(StateExp); % defines the type of experiment (0: instructed task; 1: inferred task)
InputAll(iSubject, iDate).FixationAlpha = double(FixationAlpha); % not used variable
InputAll(iSubject, iDate).tR = double(tR); % reaction time
InputAll(iSubject, iDate).Cue = double(Cue); % actual rule of environment at each trial
InputAll(iSubject, iDate).PrAn = double(PrAn); % action (-1: anti, +1: pro)
InputAll(iSubject, iDate).tDev = double(round(tDev)); % sample intervals
InputAll(iSubject, iDate).RuleChoice = double(RuleChoice); % rule choice by subject
InputAll(iSubject, iDate).TF = double(TF) .* (double(Cue) == double(RuleChoice)); % Feedback (accuracy of both action and rule)
InputAll(iSubject, iDate).Self_TF = double(TF); 
InputAll(iSubject, iDate).DLocation = double(DLocation); % location of saccade
InputAll(iSubject, iDate).DevValues = [-320, -160, -80, -40, 0, 40, 80, 160, 320]; % sample interval ranges
InputAll(iSubject, iDate).tdMean = 850; % mean of dist
InputAll(iSubject, iDate).Monkey = Subject.Names{iSubject}; % name of subject

if exist('Rule_RT') %  % reaction time of rule response 
    InputAll(iSubject, iDate).Rule_RT = double(Rule_RT); % not used here
else
    InputAll(iSubject, iDate).Rule_RT = zeros(1, length(TF)); % not used here
end

if exist('BStimStatus') % Stimulation sessions
    InputAll(iSubject, iDate).BStimStatus = BStimStatus; % not used here
else
    InputAll(iSubject, iDate).BStimStatus = zeros(1, length(TF)); % not used here
end
if exist('CostEnvironment') % Costly environment sessions
    InputAll(iSubject, iDate).CostEnvironment = CostEnvironment; % not used here
else
    InputAll(iSubject, iDate).CostEnvironment = zeros(1, length(TF)); % not used here
end




%% Fit psychometric function.

iExp = 1; % instructed rule task (rule is cued)

% Step #1: find trials that are in (instructed rule task) and make an array of those trials
Index_of_trials_forFit = find( (InputAll(iSubject, iDate).StateExp==0) );
SelectedInput.tDev = InputAll(iSubject, iDate).tDev(Index_of_trials_forFit);
SelectedInput.Cue = InputAll(iSubject, iDate).Cue(Index_of_trials_forFit);
SelectedInput.RuleChoice = InputAll(iSubject, iDate).RuleChoice(Index_of_trials_forFit);
SelectedInput.TF = InputAll(iSubject, iDate).TF(Index_of_trials_forFit);
SelectedInput.PrAn = InputAll(iSubject, iDate).PrAn(Index_of_trials_forFit);
SelectedInput.tdMean = InputAll(iSubject, iDate).tdMean;
SelectedInput.DevValues = InputAll(iSubject, iDate).DevValues;

% Step #2: fitting the psychometric model (Subjective rule)
SubjObj_flag = 'Subj'; % this fits based on the subject rule tags (see paper for more info)
[EstimatedParameters, syntheticInput, parameter_strings] = ModelParameterOptimization(SelectedInput, SubjObj_flag);
model.psychometric.subj.p_anti(1, :, :) = syntheticInput.p_anti; % fitted curve
model.psychometric.subj.x_axis = syntheticInput.x_axis; % x axis of fit for psychometric function
model.psychometric.subj.EstimatedParameters(1, :) = EstimatedParameters; % fitted parameters
model.psychometric.subj.parameter_strings = parameter_strings; % name of fitted parameters
model.psychometric.subj.SubjObj_flag = SubjObj_flag; % a tags shows subjective or objective fit

if test_plot_enable == 1 % plot the subjective psychometric function for the instructed task
   temp(:,:) = model.psychometric.subj.p_anti(1, :, :);
   x_axis = model.psychometric.subj.x_axis;
   figure; plot(x_axis, temp'); lH = legend('RuleA', 'RuleB'); set(lH, 'Box', 'off'); xlabel('td'); ylabel('Pr(Anti | rule, td)'); title('Subjective psychometric  in instructed rule task');
    
end

% Step #3: fitting the psychometric model (Objective rule) and make an array of those trials
SubjObj_flag = 'Obj'; % this fits based on the objective rule tags (see paper for more info)
[EstimatedParameters, syntheticInput, parameter_strings] = ModelParameterOptimization(SelectedInput, SubjObj_flag);
model.psychometric.obj.p_anti(1, :, :) = syntheticInput.p_anti; % fitted curve
model.psychometric.obj.x_axis = syntheticInput.x_axis;
model.psychometric.obj.EstimatedParameters(1, :) = EstimatedParameters;
model.psychometric.obj.parameter_strings = parameter_strings;
model.psychometric.obj.SubjObj_flag = SubjObj_flag;

if test_plot_enable == 1 % plot the objective psychometric function for the instructed rule task
   temp(:,:) = model.psychometric.subj.p_anti(1, :, :);
   x_axis = model.psychometric.obj.x_axis;
   figure; plot(x_axis, temp'); lH = legend('RuleA', 'RuleB'); set(lH, 'Box', 'off'); xlabel('td'); ylabel('Pr(Anti | rule, td)'); title('Objective psychometric  in instructed rule task');
end


% Step #4: find trials that are in (inferred rule task)
Index_of_trials_forFit = find( (InputAll(iSubject, iDate).StateExp==1) );
% make an array of those trials
SelectedInput.tDev = InputAll(iSubject, iDate).tDev(Index_of_trials_forFit);
SelectedInput.Cue = InputAll(iSubject, iDate).Cue(Index_of_trials_forFit);
SelectedInput.RuleChoice = InputAll(iSubject, iDate).RuleChoice(Index_of_trials_forFit);
SelectedInput.TF = InputAll(iSubject, iDate).TF(Index_of_trials_forFit);
SelectedInput.PrAn = InputAll(iSubject, iDate).PrAn(Index_of_trials_forFit);
SelectedInput.tdMean = InputAll(iSubject, iDate).tdMean;
SelectedInput.DevValues = InputAll(iSubject, iDate).DevValues;

% Step #5: fitting the psychometric model (Subjective rule)
SubjObj_flag = 'Subj'; % this fits based on the subject rule tags (see paper for more info)
[EstimatedParameters, syntheticInput, parameter_strings] = ModelParameterOptimization(SelectedInput, SubjObj_flag);
model.psychometric.subj.p_anti(2, :, :) = syntheticInput.p_anti; % fitted curve
model.psychometric.subj.x_axis = syntheticInput.x_axis; % x axis of fit for psychometric function
model.psychometric.subj.EstimatedParameters(2, :) = EstimatedParameters; % fitted parameters
model.psychometric.subj.parameter_strings = parameter_strings; % name of fitted parameters
model.psychometric.subj.SubjObj_flag = SubjObj_flag; % a tags shows subjective or objective fit

if test_plot_enable == 1 % plot the subjective psychometric function for the inferred task
   temp(:,:) = model.psychometric.subj.p_anti(2, :, :);
   x_axis = model.psychometric.obj.x_axis;
   figure; plot(x_axis, temp'); lH = legend('RuleA', 'RuleB'); set(lH, 'Box', 'off'); xlabel('td'); ylabel('Pr(Anti | rule, td)'); title('Objective psychometric  in inferred rule task');
end


%note: psychometric model has two options for scalar and non-scalar model.
%and it can be selected in the psychometric and optimizer code.

%% fit the belief update model (this is after we fitted the psychometric functions):

% step #6: Sort the trials based on the coditions (1-B error,  2B error) in the inferred task
temp_parameters_as_input(:) = model.psychometric.subj.EstimatedParameters(2, :); % Subjective psychometric parameters. Index can be 1 (instructed task) or 2 (inferred task)
temp_parameter_strings = model.psychometric.subj.parameter_strings; % name of fitted parameters
temp_SubjObj_flag = model.psychometric.subj.SubjObj_flag; % The tag shows subjective or objective fit

    % prepare data format to feed the modelling sub-function
        % extract 1-B error trials (1-B Error and 2-B Reward)
            conditions_1Back_index = find( (InputAll(iSubject, iDate).StateExp(1:end-2) == 1) .* (InputAll(iSubject, iDate).TF(1:end-2)==1) .* (InputAll(iSubject, iDate).TF(2:end-1)==0) ) +1; % RW, ER  => {1B-Er}
            conditions_2Back_index = find( (InputAll(iSubject, iDate).StateExp(1:end-3) == 1) .* (InputAll(iSubject, iDate).TF(1:end-3)==1) .* (InputAll(iSubject, iDate).TF(2:end-2)==0) .* (InputAll(iSubject, iDate).TF(3:end-1)==0) ) +2; % RW, ER, ER  => {2B-Er}
            
            % making a dataset for the modelling sub-function:
            for iTrial =1: length(conditions_1Back_index)
                switchInput(iTrial).T = 1; % number of 1Back errors
                switchInput(iTrial).tDev = InputAll(iSubject, iDate).tDev(conditions_1Back_index(iTrial)); % stimulus
                switchInput(iTrial).PrAn = InputAll(iSubject, iDate).PrAn(conditions_1Back_index(iTrial)); % Action
                switchInput(iTrial).RuleChoice = InputAll(iSubject, iDate).RuleChoice(conditions_1Back_index(iTrial)); % rule choice
                switchInput(iTrial).Cue = InputAll(iSubject, iDate).Cue(conditions_1Back_index(iTrial)); % actual rule
                switchInput(iTrial).TF = InputAll(iSubject, iDate).TF(conditions_1Back_index(iTrial)); % feedback
                switchInput(iTrial).SW = InputAll(iSubject, iDate).RuleChoice(conditions_1Back_index(iTrial)) ~= InputAll(iSubject, iDate).RuleChoice(conditions_1Back_index(iTrial)+1); % switch/notswitch
                switchInput(iTrial).DevValues = InputAll(iSubject, iDate).DevValues; % array of stimulus (samples)
                switchInput(iTrial).tdMean = InputAll(iSubject, iDate).tdMean; % mean of sample distribution
            end
            
            for iTrial = 1: length(conditions_2Back_index)
                switchInput(length(conditions_1Back_index)+iTrial).T = 2; % number of 1Back errors
                switchInput(length(conditions_1Back_index)+iTrial).tDev = [InputAll(iSubject, iDate).tDev(conditions_2Back_index(iTrial) -1), InputAll(iSubject, iDate).tDev(conditions_2Back_index(iTrial))];
                switchInput(length(conditions_1Back_index)+iTrial).PrAn = [InputAll(iSubject, iDate).PrAn(conditions_2Back_index(iTrial) -1), InputAll(iSubject, iDate).PrAn(conditions_2Back_index(iTrial))];
                switchInput(length(conditions_1Back_index)+iTrial).RuleChoice = [InputAll(iSubject, iDate).RuleChoice(conditions_2Back_index(iTrial) -1), InputAll(iSubject, iDate).RuleChoice(conditions_2Back_index(iTrial))];
                switchInput(length(conditions_1Back_index)+iTrial).Cue = [InputAll(iSubject, iDate).Cue(conditions_2Back_index(iTrial) -1), InputAll(iSubject, iDate).Cue(conditions_2Back_index(iTrial))];
                switchInput(length(conditions_1Back_index)+iTrial).TF = [InputAll(iSubject, iDate).TF(conditions_2Back_index(iTrial) -1), InputAll(iSubject, iDate).TF(conditions_2Back_index(iTrial))];
                
                switchInput(length(conditions_1Back_index)+iTrial).SW = InputAll(iSubject, iDate).RuleChoice(conditions_2Back_index(iTrial)) ~= InputAll(iSubject, iDate).RuleChoice(conditions_2Back_index(iTrial)+1);
                switchInput(length(conditions_1Back_index)+iTrial).DevValues = InputAll(iSubject, iDate).DevValues; 
                switchInput(length(conditions_1Back_index)+iTrial).tdMean = InputAll(iSubject, iDate).tdMean; 
            end

     % Fit the switch model
            Set_param_variable = [];
            [~, ~, estimate_parameters, ~, ~, MachineSimulation] = SwitchParameterOptimization(switchInput, temp_parameters_as_input, temp_parameter_strings, temp_SubjObj_flag, InputAll(iSubject, iDate), Set_param_variable);       % pr_switch_model(iRuleChoice, itd, T_numOfBackError )
            % model parameters:
            % 1- estimate_parameters.sigma_switch_estimated (sigma parameter of the switch model) (Free parameter for optimization)
            % 2- estimate_parameters.pam3_estimated (alpha of the switch model) (Free parameter for optimization)
            % hazard parameter is fixed.
            % consecutive errors (for example. for T=2 => sigma2= sigma * sqrt(2), ...)
                        
      % measuring switch probability of subject and also the machineSimulation
            % for subject:
                for itd = 1: length(InputAll(iSubject, iDate).DevValues) % for each sample interval:
                    % for 1-Back:
                    conditions_1Back_index = find( (InputAll(iSubject, iDate).StateExp(1:end-2) == 1) .* (InputAll(iSubject, iDate).TF(1:end-2)==1) .* (InputAll(iSubject, iDate).TF(2:end-1)==0) .* (InputAll(iSubject, iDate).tDev(2:end-1)==InputAll(iSubject, iDate).DevValues(itd)) ) +1; % RW, ER  => {1B-Er}
                    SW_v1 = InputAll(iSubject, iDate).RuleChoice(conditions_1Back_index) ~= InputAll(iSubject, iDate).RuleChoice(conditions_1Back_index+1); 
                    iT = 1; Pr_Sw(iT, itd) = sum(SW_v1) / length(SW_v1);
                    % for 2-Back Er (cosecutive errors with no 1-B switch):
                    conditions_2Back_index = find( (InputAll(iSubject, iDate).StateExp(1:end-3) == 1) .* (InputAll(iSubject, iDate).TF(1:end-3)==1) .* (InputAll(iSubject, iDate).TF(2:end-2)==0) .* (InputAll(iSubject, iDate).TF(3:end-1)==0) .* (InputAll(iSubject, iDate).RuleChoice(2:end-2) == InputAll(iSubject, iDate).RuleChoice(3:end-1)) .* (InputAll(iSubject, iDate).tDev(3:end-1)==InputAll(iSubject, iDate).DevValues(itd)) ) +2; % RW, ER, ER  => {2B-Er}
                    SW_v2 = InputAll(iSubject, iDate).RuleChoice(conditions_2Back_index) ~= InputAll(iSubject, iDate).RuleChoice(conditions_2Back_index+1); 
                    iT = 2; Pr_Sw(iT, itd) = sum(SW_v2) / length(SW_v2);
                end
            % for machine simulation:
                for itd = 1: length(InputAll(iSubject, iDate).DevValues) % for each sample interval:
                    % for 1-Back:
                    conditions_1Back_index = find( (MachineSimulation{1}.StateExp(1:end-2) == 1) .* (MachineSimulation{1}.TF(1:end-2)==1) .* (MachineSimulation{1}.TF(2:end-1)==0) .* (MachineSimulation{1}.tDev(2:end-1)==InputAll(iSubject, iDate).DevValues(itd)) ) +1; % RW, ER  => {1B-Er}
                    SW_v1 = MachineSimulation{1}.RuleChoice(conditions_1Back_index) ~= MachineSimulation{1}.RuleChoice(conditions_1Back_index+1); 
                    iT = 1; Pr_Sw_m(iT, itd) = sum(SW_v1) / length(SW_v1);
                    % for 2-Back Er (cosecutive errors with no 1-B switch):
                    conditions_2Back_index = find( (MachineSimulation{1}.StateExp(1:end-3) == 1) .* (MachineSimulation{1}.TF(1:end-3)==1) .* (MachineSimulation{1}.TF(2:end-2)==0) .* (MachineSimulation{1}.TF(3:end-1)==0) .* (MachineSimulation{1}.RuleChoice(2:end-2) == MachineSimulation{1}.RuleChoice(3:end-1)) .* (MachineSimulation{1}.tDev(3:end-1)==InputAll(iSubject, iDate).DevValues(itd)) ) +2; % RW, ER, ER  => {2B-Er}
                    SW_v2 = MachineSimulation{1}.RuleChoice(conditions_2Back_index) ~= MachineSimulation{1}.RuleChoice(conditions_2Back_index+1); 
                    iT = 2; Pr_Sw_m(iT, itd) = sum(SW_v2) / length(SW_v2);
                end
                
            figure; hold on;
            iT = 1; plot(InputAll(iSubject, iDate).DevValues, Pr_Sw(iT, :), 'o', 'Color', [1,0,0], 'LineWidth', 1.5); % subject , 1B error
            iT = 2; plot(InputAll(iSubject, iDate).DevValues, Pr_Sw(iT, :), 'o', 'Color', [0.5,0,0], 'LineWidth', 1.5); % subject , 2B error
            
            iT = 1; plot(InputAll(iSubject, iDate).DevValues, Pr_Sw_m(iT, :), '-', 'Color', [1,0,0],  'LineWidth', 1.5); % machine simulation (model) , 1B error
            iT = 2; plot(InputAll(iSubject, iDate).DevValues, Pr_Sw_m(iT, :), '-', 'Color', [0.5,0,0],  'LineWidth', 1.5); % machine simulation (model) , 2B error
            ylim([0,1]);
            ylabel('Pr(SW)'); xlabel('td');
            legend('Subject, 1B-Er', 'Subject, 2B-Er', 'Model, 1B-Er', 'Model, 2B-Er');
            
            
            
    
