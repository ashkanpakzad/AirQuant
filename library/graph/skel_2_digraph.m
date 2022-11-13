function [digraphout, glinkout, gnode] = skel_2_digraph(skel, method)
    % Generate digraph from skeleton.
    %
    % Generate digraph from skeleton using skel2graph library.
    % It is necessary for AirQuant to use a digraph as it dictates the direction
    % that splines are interpolated for tubes. The direction is
    % set by designating an origin node from which all directions in
    % the digraph are oriented away. :attr:`method` `topnode` selects the
    % highest z-point. The nodes are ordered by a BFS from the origin node.
    %
    % See https://github.com/phi-max/skel2graph3d-matlab for more
    % details on skel2graph.
    %
    %
    % Args:
    %   skel: skeleton to turn into digraph 
    %   method(string): *OPTIONAL* method to choose originating direction.
    %       default uses `topnode` method by chossing the most
    %       superior node.
    %
    % Return:
    %   3 items
    %
    %   - **digraph** : converted digraph object
    %   - **glink** (:class:`struct`): positional information for each branch relating to
    %     skel.
    %   - **gnode** (:class:`struct`): positional information for each node relating to
    %     skel.
    %
    % Example:
    %
    %   >>> [digraphout, glink, gnode] = skel_2_digraph(skel, 'topnode')
    %


    if nargin < 2
        method = 'topnode';
    end

    [gadj,gnode,glink] = Skel2Graph3D(skel,1);
    % choose originating node using chosen method


    % Create digraph with edges in both directions, loop through
    % branches and remove opposing direction to originating node.
    G = digraph(gadj);

    assert(length(glink) == height(G.Edges)/2, 'Skeleton appears to contain loops. This is not supported.')

    bins = conncomp(G);
    
    if isnumeric(method)
        assert(length(unique(bins)) < 2, ['Can only be used with one ' ...
            'connected component.'])
        assert(length(method)==3,'Must be of length 3')
        % construct list of nodes with origins into array
        node_coords = [[gnode.comy]; [gnode.comx]; [gnode.comz]]';
        originnode = dsearchn(node_coords,method);

    elseif strcmp(method,'topnode')
            % use most superiour node as origin for each subgraph
            originnode = zeros(max(bins),1);
            for ii = 1:max(bins)
                binbool = (bins==ii);
                nodelist = 1:numnodes(G);
                binidx = nodelist(binbool);
                gnodeii = gnode(binidx);
                [~, binorigin] = max([gnodeii.comz]);
                originnode(ii) = binidx(binorigin);
            end
    else
        error('Invalid method')
    end


    % BF search from carina node to distal for each disconnected graph.
    allnode_discovery = cell(max(bins),1);
    for ii = 1:max(bins)
        allnode_discovery{ii} = bfsearch(G,originnode(ii));
    end
    node_discovery = cell2mat(allnode_discovery);
    %%% reorder nodes by bfs
    G = reordernodes(G, node_discovery);
    % reorder in gnodes and glinks
    gnode = gnode(node_discovery);
    for iinode = 1:length(gnode)
        conns = gnode(iinode).conn;
        for iiconn = 1:length(conns)
            gnode(iinode).conn(iiconn) = find(node_discovery == conns(iiconn));
        end
    end
    for iilink = 1:length(glink)
        glink(iilink).n1 = find(node_discovery == glink(iilink).n1);
        glink(iilink).n2 = find(node_discovery == glink(iilink).n2);
    end

    % organise into outward facing digraph
    % half of the edges will be removed
    removal = zeros(height(G.Edges)/2 , 1);
    j = 1;
    for i = 1:height(G.Edges)
        % edges should be connected by smaller to larger node.
        if (G.Edges.EndNodes(i,1) - G.Edges.EndNodes(i,2)) > 0
            % record if not central-->distal
            removal(j) = i;
            j = j + 1;
        end
    end

    % remove recorded edges at end in opposite direction.
    removal(removal == 0) = [];
    G = rmedge(G,removal);

    % Need to ensure directions in glink are also in the correct
    % direction. Swap over if not.
    for i = 1:length(glink)
        glink_nodes = [glink(i).n1, glink(i).n2];
        % if not central-->distal
        if ~ismember(glink_nodes, G.Edges.EndNodes, 'rows')
            % Swap round node assignments
            glink(i).n1 = glink_nodes(2);
            glink(i).n2 = glink_nodes(1);
            % Swap round path order
            glink(i).point = fliplr(glink(i).point);
        end
    end    

    % Copy over glink to create digraph with same edge order
    edges = [[glink.n1]', [glink.n2]'];
    weights = zeros(length(glink),1);
    for i = 1:length(glink)
        weights(i) = length(glink(i).point);
    end
    labels = [1:length(glink)]';
    Edgetable = table(edges,weights,labels,'VariableNames',{'EndNodes', 'Weight', 'Label'});
    
    digraphout = digraph(Edgetable);
    % add node properties from Gnode
    digraphout.Nodes.comx(:) = [gnode(:).comx];
    digraphout.Nodes.comy(:) = [gnode(:).comy];
    digraphout.Nodes.comz(:) = [gnode(:).comz];
    digraphout.Nodes.ep(:) = [gnode(:).ep];
    digraphout.Nodes.label(:) = [1:length(gnode)]';

    % BFS per graph in digraph
    % first convert digraph to graph to get CCs
    bins = conncomp(digraphout,'Type','weak','OutputForm','cell');
    
    % order bins by descending size of each CC
    [~,I] = sort(cellfun(@length,bins));
    bins = bins(I);
    
    % in that order get edge BFS of each CC
    E = [];
    for ii = 1:length(bins)
        nodes = bins{ii};
        % find root, defined by min node idx
        root_node = min(nodes);
        % find sub graph (CC) BFS
        [~,sub_E] = bfsearch(digraphout,root_node,'edgetonew');
        E = [E; sub_E];
    end
    
    assert(length(E)==length(glink),'Number of edges in new edge BFS not expected.')

    % convert E indicies to glink indicies
    E_glink = digraphout.Edges.Label(E);
    
    % reorder glink
    glinkout = glink(E_glink);
end
