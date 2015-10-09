function stacks = stackGenerator(nAngles,nPlies,scenario)
  % Generate stacks of plies according to different scenarios
  % Input:
    % Number of equally spaced angles in stack: nAngles
    % Number of plies in laminate: nPlies
    % Considered scenario: scenario = 1,2 or 3
  % Output:
    % Structure containing stack configurations: stacks

  stacks = zeros(nAngles,nPlies,2);

  if scenario==1
    % 0 to 45 degrees
    angles = [0:pi/(4*(nAngles-1)):pi/4]';
    % Stacks of plies in same direction
    stacks(:,:,1) = repmat(angles,[1,nPlies]);
    % Stacks of \pm angle plies
    stacks(:,1,2)=angles;
    stacks(:,:,2) = repmat(stacks(:,1,2),[1,nPlies]);
    stacks(:,2:2:end,2) = -1.*stacks(:,2:2:end,2);

  elseif scenario==2
    maxAngle = pi/(4*nAngles):pi/(4*nAngles):pi/4;
    for i = 1:nAngles
      angles = [maxAngle(i)/((nPlies-1)/2):maxAngle(i)/((nPlies-1)/2):maxAngle(i)];
      stacks(i,:,1) = [fliplr(angles), 0, angles];
      stacks(i,:,2) = [fliplr(-1.*angles), 0, angles];
    end

  elseif scenario==3
      maxAngle = pi/(4*nAngles):pi/(4*nAngles):pi/4;
      for i = 1:nAngles
        angles = [maxAngle(i)/((nPlies-1)/2):maxAngle(i)/((nPlies-1)/2):maxAngle(i)];
        angles2 = [maxAngle(i)/((nPlies-1)/4):maxAngle(i)/((nPlies-1)/4):maxAngle(i)];
        angles2 = repelem(angles2,2);
        angles2(2:2:end) = -angles2(2:2:end);
        stacks(i,:,1) = [angles, 0, fliplr(angles)];
        stacks(i,:,2) = [angles, 0, fliplr(-1.*angles)];
        stacks(i,:,3) = [angles2, 0, fliplr(angles2)];
        stacks(i,:,4) = [fliplr(angles2), 0, angles2];
      end
  else
    error('Scenario is not yet defined!')
  end
