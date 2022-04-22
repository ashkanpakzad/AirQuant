%% Main function to generate tests
function tests = TestSkel2Digraph()
tests = functiontests(localfunctions);
end

%% fixtures  
function setupOnce(testCase)
    AirQuantAddPath();
    skel_name = 'FIXskel.nii.gz';
    testCase.TestData.skel = logical(niftiread(skel_name));
    [Gadj,~,Glink] = Skel2Graph3D(testCase.TestData.skel,0);
    G = graph(Gadj);
    testCase.TestData.Glink = Glink;
    testCase.TestData.G = G;
end

function testNumBranches(testCase)
    % same number of branches as skel2graph
    [digraphout, ~] = Skel2Digraph(testCase.TestData.skel, 'topnode');
    actual = height(digraphout.Edges);
    expected = height(testCase.TestData.G.Edges);
    verifyEqual(testCase,actual,expected)
end

function testTopNodeNoIn(testCase)
    % topnode has no in edges
    [digraphout, ~] = Skel2Digraph(testCase.TestData.skel, 'topnode');
    actual = inedges(digraphout,1);
    verifyEmpty(testCase,actual)
end

function teardownOnce(testCase)
end
