%% Main function to generate tests
function tests = Test_workflow()
tests = functiontests(localfunctions);
end

%% fixtures  
function setupOnce(testCase)
    AQpath = AirQuantAddPath;

    % download test data if necessary
    casename='chestct';
    AQdownload_data(casename);

    % set names
    data_dir = fullfile(AQpath,'data','airquant');
    results_dir = fullfile(AQpath,'data','results');

    % clear results dir if anything stored
    if exist(fullfile(results_dir,casename),'dir')
        rmdir(fullfile(results_dir,casename),'s')
    end

    % get file names
    filesuffix = '.nii.gz';
    sourcef = fullfile(data_dir,[casename,'_source',filesuffix]);
    segf = fullfile(data_dir,[casename,'_airway',filesuffix]);
    skelf = fullfile(data_dir,[casename,'_airway_PTKskel',filesuffix]);

    % save vars
    testCase.TestData = struct('casename',casename,'sourcef',sourcef,'segf',segf,'skelf',skelf,'results_dir',results_dir);
end

function test_runcases(testCase)
    vars=testCase.TestData;
    wf_clinicalairways_fwhmesl(vars.casename, vars.sourcef, vars.segf, vars.skelf, vars.results_dir);
    % check AQnet saves
    verifyTrue(testCase,exist(fullfile(vars.results_dir,vars.casename,[vars.casename,'_AQnet.mat']),'file')>0)
end

function teardownOnce(testCase)
end
