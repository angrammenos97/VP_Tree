clear
n = 10; dim = 2; k = 3; % works for any d
X = rand(n,dim);
% assume vantage point is the last one
% get squares of distances from it
d = sqrt( sum((X(1:n-1,:) - X(n,:)).^2,2) );
% and find the median distance
medianDistance = median(d);
% plot them to confirm
clf; hold on; axis equal
plot(X(d <= medianDistance,1), X(d <= medianDistance,2), 'r.')
plot(X(d > medianDistance,1), X(d > medianDistance,2), 'b.')
plot(X(n,1), X(n,2), 'bo') % vantage point
% create the search tree
tree = vpTree(X);
% search the tree
q = rand(1,dim); % random query point
plot(q(1,1), q(1,2), 'go')
tau = sqrt(2); % desired search radius
list = repmat([0 tau], k, 1); % list with all neighbors
list = searchNb(tree, list, tau, q);

function T = vpTree(X)
    % function T = vpTree(X)
    % computes the vantage point tree structure with
    % T.vp : the vantage point
    % T.md : the median distance of the vantage point to the other points
    % T.idx : the index of the vantage point in the original set
    % T.inner and T.outer : vantage-point subtrees
    % of a set of n points in d dimensions with coordinates X(1:n,1:d)
    %
    T = vpt(X, 1:size(X,1));
    
    function T = vpt(X, idx)
        n = size(X,1); % number of points
        if n == 0
            T = [];
        else
            T.vp = X(n,:);
            T.idx = idx(n);
            d = sqrt( sum((X(1:n-1,:) - X(n,:)).^2,2) );
            medianDistance = median(d);
            T.md = medianDistance;
            % split and recurse
            inner = d <= medianDistance;
            T.inner = vpt(X( inner,:), idx( inner));
            T.outer = vpt(X(~inner,:), idx(~inner));
        end
    end
end

function list = searchNb(tree, list, tau, q)
    if (isempty(tree) )
        return;
    else
        dist = sqrt( sum((tree.vp(1,:) - q(1,:)).^2,2) );        
        if (dist < tau)
            tauId = find(list(:,2) == tau);
            list(tauId(1),:) = [tree.idx, dist];
            M = max(list);
            tau = M(1,2);
        end
        if (dist < tree.md)
            list = searchNb(tree.inner, list, tau,q);
            if (tree.md < (dist + tau))
                list = searchNb(tree.outer, list, tau, q);
            end
        else
            list = searchNb(tree.outer, list, tau,q);
            if (tree.md > abs(dist - tau))
                list = searchNb(tree.inner, list, tau, q);
            end
        end       
    end    
end