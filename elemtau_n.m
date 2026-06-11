function [ tau , normal ] = elemtau_n( node , elem , elem2edgeSign ,elem2edge ,edge )
NT = size(elem,1);
tau = zeros(NT,2,3);
normal = tau;
tau(:,:,1) = (node(edge(elem2edge(:,1),2),:) - node(edge(elem2edge(:,1),1),:)).*elem2edgeSign(:,1);
tau(:,:,2) = (node(edge(elem2edge(:,2),2),:) - node(edge(elem2edge(:,2),1),:)).*elem2edgeSign(:,2);
tau(:,:,3) = (node(edge(elem2edge(:,3),2),:) - node(edge(elem2edge(:,3),1),:)).*elem2edgeSign(:,3);
tau(:,:,1) = tau(:,:,1)./sqrt(tau(:,1,1).^2+tau(:,2,1).^2);
tau(:,:,2) = tau(:,:,2)./sqrt(tau(:,1,2).^2+tau(:,2,2).^2);
tau(:,:,3) = tau(:,:,3)./sqrt(tau(:,1,3).^2+tau(:,2,3).^2);
normal(:,:,1) = [tau(:,2,1),-tau(:,1,1)];
normal(:,:,2) = [tau(:,2,2),-tau(:,1,2)];
normal(:,:,3) = [tau(:,2,3),-tau(:,1,3)];
end

