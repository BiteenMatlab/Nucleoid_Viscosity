<<<<<<< HEAD
function Result = BD_Simu_batchGen(Num_run,StepSize,result_folder,scale_bar,Ro,Lo,kb_nuc_spMap,anchor_spMap,RndStep,ConfMethod,num_step,prog_bar)
% BD_Simu_batchGen performs a batch of Brownian dynamics simulations in 3D.
%
% INPUTS:
%   Num_run       - Number of simulation runs.
%   StepSize      - Step size in nm.
%   result_folder - Folder (string or char) in which the simulation results
%                   will be stored. When empty, nothing get stored. 
%   scale_bar     - Scale factor (nonzero numeric scalar) to adjust StepSize.
%   Ro            - Parameter for BD_Simu_3D.
%   Lo            - Parameter for BD_Simu_3D.
%   kb_nuc_spMap  - Parameter for BD_Simu_3D.
%   anchor_spMap  - Parameter for BD_Simu_3D.
%   RndStep       - (Optional) 'on' or 'off', default = 'on'.
%   ConfMethod    - (Optional) confinement method string, default = 'Collision'.
%   num_step      - (Optional) number of steps for each trajectory, default = 1e5.
%   prog_bar      - (Optional) show progress bar or not. 
%
% Author: Xiaofeng Dai
% Date: 05/22/2025

if nargin < 8
    error('Not enough input arguments. At least 8 inputs are required.');
end
if ~isscalar(Num_run) || ~isnumeric(Num_run) || (Num_run < 1)
    error('Num_run must be a positive numeric scalar.');
end
if ~isscalar(StepSize) || ~isnumeric(StepSize) || (StepSize <= 0)
    error('StepSize must be a positive numeric scalar.');
end
if ~(ischar(result_folder) || isstring(result_folder))   &&  ~isempty(result_folder)
    error('result_folder must be a string or a character array.');
end

% Set default optional parameters if necessary  
if nargin == 8
    RndStep = 'on';
    ConfMethod = 'Collision';
    num_step = 1e5;
    prog_bar = 'off';
end
if nargin == 9
    ConfMethod = 'Collision';
    num_step = 1e5;
    prog_bar = 'off';
end
if nargin == 10
    num_step = 1e5;
    prog_bar = 'off';
end
if nargin == 11
    prog_bar = 'off';
end

if ~isempty(result_folder)
    % Create main folder to store Brownian dynamics simulation results
    BDSimu_folder = [result_folder,filesep,'BDSimu Results'];
    if ~exist(BDSimu_folder, 'dir')
    mkdir(BDSimu_folder)
    end
    
    % Create subfolder for storing results with the given step size
    StepStore_folder = strcat(BDSimu_folder,filesep,num2str(StepSize),'nm');
    if ~exist(StepStore_folder, 'dir')
    mkdir(StepStore_folder)
    end
end

% Preallocate a cell array to hold simulation results for each run
Result = cell(Num_run,1);

%--------------------------------------------------------------------------
% Set up a progress bar for the simulation runs in the parallel loop
useProgBar = strcmpi(prog_bar,'on');
% placeholder so it always exists (for parafor)
progressQueue = [];      
nCompleted    = 0;       
hWaitBar      = [];      
if useProgBar
    hWaitBar = waitbar(0, strcat('BD simulations, step size = ',string(StepSize),' nm ...'));
    nCompleted = 0;
    progressQueue = parallel.pool.DataQueue;
    afterEach(progressQueue, @updateWaitBar);
end
%--------------------------------------------------------------------------

% Run the simulations in parallel using a parfor loop
parfor i = 1:Num_run
    try
        % Run a single Brownian dynamics simulation for 1e5 steps.
        Traj = BD_Simu_3D(num_step, StepSize/scale_bar, Ro, Lo, ...
            kb_nuc_spMap, anchor_spMap, ...
            'RandomStep', RndStep, 'Confinement', ConfMethod, ...
            'PlotTraj', 'off');
        Result{i} = Traj;
    catch ME
        % If an error occurs, note the run (displaying i) and issue a warning with the error message.
        warning(['Empty in ', num2str(i), ' run: ', ME.message]);
    end
    if useProgBar
        % Send a signal to update the progress bar regardless of success/failure.
        send(progressQueue, []);
    end
end

if useProgBar && isvalid(hWaitBar)
    % Once the simulation runs are done, close the progress bar.
    delete(hWaitBar)
end

clear i Traj ind

%--------------------------------------------------------------------------
if ~isempty(result_folder)
    % Save each simulation result to a .mat file.
    padLen = floor(log10(Num_run)) + 1;
    for i = 1:Num_run
        store = Result{i};
        Number = sprintf(['%0', num2str(padLen), 'd'], i);
        filename1 = strcat('Simu_result_', num2str(StepSize), 'nm_', Number, '.mat');
        file1 = strcat(StepStore_folder, filesep, filename1);
        save(file1, "store")
    end
end
%--------------------------------------------------------------------------

clear Num_run store index Number
clear filename1 file1 i
clear BDSimu_folder StepStore_folder

    % Nested function: Update the waitbar progress based on the number of completed iterations.
    function updateWaitBar(~)
        nCompleted = nCompleted + 1;
        waitbar(nCompleted/Num_run, hWaitBar);
    end

=======
function Result = BD_Simu_batchGen(Num_run,StepSize,result_folder,scale_bar,Ro,Lo,kb_nuc_spMap,anchor_spMap,RndStep,ConfMethod,num_step,prog_bar)
% BD_Simu_batchGen performs a batch of Brownian dynamics simulations in 3D.
%
% INPUTS:
%   Num_run       - Number of simulation runs.
%   StepSize      - Step size in nm.
%   result_folder - Folder (string or char) in which the simulation results
%                   will be stored. When empty, nothing get stored. 
%   scale_bar     - Scale factor (nonzero numeric scalar) to adjust StepSize.
%   Ro            - Parameter for BD_Simu_3D.
%   Lo            - Parameter for BD_Simu_3D.
%   kb_nuc_spMap  - Parameter for BD_Simu_3D.
%   anchor_spMap  - Parameter for BD_Simu_3D.
%   RndStep       - (Optional) 'on' or 'off', default = 'on'.
%   ConfMethod    - (Optional) confinement method string, default = 'Collision'.
%   num_step      - (Optional) number of steps for each trajectory, default = 1e5.
%   prog_bar      - (Optional) show progress bar or not. 
%
% Author: Xiaofeng Dai
% Date: 05/22/2025

if nargin < 8
    error('Not enough input arguments. At least 8 inputs are required.');
end
if ~isscalar(Num_run) || ~isnumeric(Num_run) || (Num_run < 1)
    error('Num_run must be a positive numeric scalar.');
end
if ~isscalar(StepSize) || ~isnumeric(StepSize) || (StepSize <= 0)
    error('StepSize must be a positive numeric scalar.');
end
if ~(ischar(result_folder) || isstring(result_folder))   &&  ~isempty(result_folder)
    error('result_folder must be a string or a character array.');
end

% Set default optional parameters if necessary  
if nargin == 8
    RndStep = 'on';
    ConfMethod = 'Collision';
    num_step = 1e5;
    prog_bar = 'off';
end
if nargin == 9
    ConfMethod = 'Collision';
    num_step = 1e5;
    prog_bar = 'off';
end
if nargin == 10
    num_step = 1e5;
    prog_bar = 'off';
end
if nargin == 11
    prog_bar = 'off';
end

if ~isempty(result_folder)
    % Create main folder to store Brownian dynamics simulation results
    BDSimu_folder = [result_folder,filesep,'BDSimu Results'];
    if ~exist(BDSimu_folder, 'dir')
    mkdir(BDSimu_folder)
    end
    
    % Create subfolder for storing results with the given step size
    StepStore_folder = strcat(BDSimu_folder,filesep,num2str(StepSize),'nm');
    if ~exist(StepStore_folder, 'dir')
    mkdir(StepStore_folder)
    end
end

% Preallocate a cell array to hold simulation results for each run
Result = cell(Num_run,1);

%--------------------------------------------------------------------------
% Set up a progress bar for the simulation runs in the parallel loop
useProgBar = strcmpi(prog_bar,'on');
% placeholder so it always exists (for parafor)
progressQueue = [];      
nCompleted    = 0;       
hWaitBar      = [];      
if useProgBar
    hWaitBar = waitbar(0, strcat('BD simulations, step size = ',string(StepSize),' nm ...'));
    nCompleted = 0;
    progressQueue = parallel.pool.DataQueue;
    afterEach(progressQueue, @updateWaitBar);
end
%--------------------------------------------------------------------------

% Run the simulations in parallel using a parfor loop
parfor i = 1:Num_run
    try
        % Run a single Brownian dynamics simulation for 1e5 steps.
        Traj = BD_Simu_3D(num_step, StepSize/scale_bar, Ro, Lo, ...
            kb_nuc_spMap, anchor_spMap, ...
            'RandomStep', RndStep, 'Confinement', ConfMethod, ...
            'PlotTraj', 'off');
        Result{i} = Traj;
    catch ME
        % If an error occurs, note the run (displaying i) and issue a warning with the error message.
        warning(['Empty in ', num2str(i), ' run: ', ME.message]);
    end
    if useProgBar
        % Send a signal to update the progress bar regardless of success/failure.
        send(progressQueue, []);
    end
end

if useProgBar && isvalid(hWaitBar)
    % Once the simulation runs are done, close the progress bar.
    delete(hWaitBar)
end

clear i Traj ind

%--------------------------------------------------------------------------
if ~isempty(result_folder)
    % Save each simulation result to a .mat file.
    padLen = floor(log10(Num_run)) + 1;
    for i = 1:Num_run
        store = Result{i};
        Number = sprintf(['%0', num2str(padLen), 'd'], i);
        filename1 = strcat('Simu_result_', num2str(StepSize), 'nm_', Number, '.mat');
        file1 = strcat(StepStore_folder, filesep, filename1);
        save(file1, "store")
    end
end
%--------------------------------------------------------------------------

clear Num_run store index Number
clear filename1 file1 i
clear BDSimu_folder StepStore_folder

    % Nested function: Update the waitbar progress based on the number of completed iterations.
    function updateWaitBar(~)
        nCompleted = nCompleted + 1;
        waitbar(nCompleted/Num_run, hWaitBar);
    end

>>>>>>> 5e4cf7b (Initial commit)
end