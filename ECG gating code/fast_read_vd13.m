% =============================================================================
function [KSpace, varargout] = fast_read_vd13(filenameIn,varargin)
% -----------------------------------------------------------------------------
% fast_read_vd13  Read a VD13 measurement data file with the property that all 
%                  lines for a given measurement have the same number of 
%                  samples.  This differs slightly from the VB reader 
%                  where only one measurement was stored per file.  VD can
%                  store multiple measurements.
%
% KSpace = fast_read_vd13() 
%    prompts user with gui and returns kSpace cell array, each cell storing
%               the kSpace from a particular measurement.
%
% KSpace = fast_read_vd13(filename, 'option1', 'option2', 'option3', ....)
%    reads filename and returns kspace with options applied where
%    options are strings
%
% fast_read_vd13 with no arguments or return values will list options
%
% [KSpace Other] = fast_read_vd13(filename, 'ReadOther') 
%    reads filename and returns kSpace as well as Other if it exists where
% -----------------------------------------------------------------------------

   global MDH_SCANSIZE MDH_CHANSIZE MDH_PMUSIZE MAX_BYTES_PER_FREAD
   global mdhColumns mdhBitFlag

   %% Nested functions (share parent's scope)
   function set_outputs_to_zero
      KSpace{1}.Normal.data = [];
      for iArgLocal=1:nargout-1
         varargout{iArgLocal} = [];
      end
      if ~(Control.User.noGui || Control.User.silent)
         if ~isempty(Control.GUI.hWaitBar)
            close(Control.GUI.hWaitBar);
         end
      end
   end

   if (nargin == 0) && (nargout == 0)
      UserInputDefault = define_user_input_default_struct_array;
      UserInputDefault = derive_user_input_default_fieldname(UserInputDefault);
      print_available_options(UserInputDefault);
      return
   end

   % Initializing
   KSpace{1}.Normal.data = [];
   for iArg=1:nargout-1
      varargout{iArg} = [];
   end
   nOptionsIn = nargin-1;
   if nOptionsIn > 0
      optionsIn = varargin;
   else
      optionsIn = [];
   end
   set_global_mdh_parameters(filenameIn);

   Control = get_controls_from_user_input(filenameIn, optionsIn, nargin, nargout);
   if isempty(Control)
      set_outputs_to_zero;
      return
   end

   % Swap SEG and ECO indices
   % This was introduced to handle HIFU meas dat where they used PHASCOR flags
   %    to store special navigators.  Unfortunately, they were stored in a way
   %    contrary to how PHASCOR are stored in general.  To read them, add
   %    'ReadPhaseCor' and 'SwapSegAndEco' to the call parameters.
   if Control.User.swapSegAndEco
      temp = mdhColumns.scan.indexECO;
      mdhColumns.scan.indexECO = mdhColumns.scan.indexSEG;
      mdhColumns.scan.indexSEG = temp;
      clear temp
   end

   if ~(Control.User.noGui || Control.User.silent)
      Control.GUI.hWaitBar = waitbar(0, 'Extracting TextHeader', 'Name', 'Progress', 'Resize', 'on');
   end

   fid = fopen(Control.File.In.name,'r','ieee-le');

   % Confirm that we have a vd13 file
   [isVD13, MrParcRaidFileHeader, MrParcRaidFileCell] = is_vd13_file(Control, fid);
   if ~isVD13
      set_outputs_to_zero;
      return
   end

   % Loop through measurements
   for iMeas = 1:MrParcRaidFileHeader.nMeas

      if Control.User.lastMeasDatOnly && (iMeas < MrParcRaidFileHeader.nMeas)
         KSpace{iMeas}.Normal.data = [];
         continue
      end

      % Parse the text header for dimensions, scale factors, control, etc.
      [Control, MrParcRaidFileCell, Dim, CoilSelectMap, ...
       CoilSelectMeasMap, RelativeSliceNumber, WiPMemBlock] = ...
          parse_meas_text_header(Control,fid,MrParcRaidFileCell,iMeas);
      if isempty(Control)
         set_outputs_to_zero;
         return
      end

      % First pass to chunk this measurement into blocks with the same scan lengths
      if ~(MrParcRaidFileCell{iMeas}.Control.User.noGui || MrParcRaidFileCell{iMeas}.Control.User.silent)
         waitbar(0.5, Control.GUI.hWaitBar, 'Measuring number of scans');
      end

      % With readout and PMU being different sizes, 
      %    we now have to chunk in blocks of same sized
      %    lines.  That requires scanning all line lengths
      %    up front.
      [sizeAllScansInBytes, sizeChanInBytes, MdhBlock, ...
       nBlockSYNCDATA, nBlockACQEND, nBlockScan, blockNScan, ...
       nScanInMeas, dmaLengthScan, nSamplesInScan] = chunk_into_equal_line_length_blocks(fid,MrParcRaidFileCell, iMeas, Dim, Control);
      if isempty(sizeAllScansInBytes)
         set_outputs_to_zero;
         return
      end

      % Second pass, Part 1 to extract PMU
      if (nBlockSYNCDATA > 0) && MrParcRaidFileCell{iMeas}.Control.User.readPMU
         PMUOutCell{iMeas} = extract_pmu(fid,MdhBlock);
         if isempty(PMUOutCell{iMeas})
            set_outputs_to_zero;
            return
         end
      end % if we're loading PMU

      % Second pass, Part 2 to assemble dimensions and read order for readouts
      [dimNames, nScan, n, kSpaceCentreLineNo, kSpaceCentrePartitionNo, ...
         order, dimOut, dimOutProd, channelIDUniq] = ...
         scan_mdh_for_dim_and_order(fid,Dim,MrParcRaidFileCell{iMeas}.Control,sizeChanInBytes,MdhBlock,nScanInMeas,nSamplesInScan);
      if isempty(dimNames)
         set_outputs_to_zero;
         return
      end

      if ~(MrParcRaidFileCell{iMeas}.Control.User.noGui || MrParcRaidFileCell{iMeas}.Control.User.silent)
         waitbar(1.0, Control.GUI.hWaitBar, 'Output dimensions computed.');
      end

      % Third pass to read the data
      if ~(MrParcRaidFileCell{iMeas}.Control.User.noGui || MrParcRaidFileCell{iMeas}.Control.User.silent)
         waitbar(0.0, Control.GUI.hWaitBar, 'Allocating arrays');
      end

      % Allocate arrays
      progressPercent = 0.0;
      for typeCell = {'NoiseAdj' 'PhaseCor' 'RTFeedback' 'Normal'}
         typeStr = typeCell{1};
         foundStr = strcat('found',typeStr);
         if nScan.(typeStr) > 0
            Control.User.(foundStr) = true;
            if strcmp(typeStr,'PhaseCor')
               KSpace{iMeas}.(typeStr).data = complex(single(0) * zeros(n.COL, compute_product(dimOut.(typeStr)), n.SEG.(typeStr), 'single'));
               if Control.User.readTimeStamp
                  KSpace{iMeas}.(typeStr).TimeStamp.count = 0 * zeros(compute_product(dimOut.(typeStr)), n.SEG.(typeStr));
                  KSpace{iMeas}.(typeStr).TimeStamp.value = uint32(0) * zeros(compute_product(dimOut.(typeStr)), n.SEG.(typeStr), 1,'uint32');
               end
            else
               KSpace{iMeas}.(typeStr).data = complex(single(0) * zeros(n.COL, compute_product(dimOut.(typeStr)), 'single'));
               if n.SEG.(typeStr) > 1
                  KSpace{iMeas}.(typeStr).SEG = uint16(0) * zeros(1,compute_product(dimOut.(typeStr)),'uint16');
               end
               if Control.User.readTimeStamp
                  KSpace{iMeas}.(typeStr).TimeStamp.count = 0 * zeros(compute_product(dimOut.(typeStr)),1);
                  KSpace{iMeas}.(typeStr).TimeStamp.value = uint32(0) * zeros(compute_product(dimOut.(typeStr)),1,'uint32');
               end
            end
         else
            if Control.User.readPhaseCor & strcmp(typeStr,'PhaseCor')
               KSpace{iMeas}.(typeStr).data = [];
            end
            if Control.User.readNoiseAdj & strcmp(typeStr,'NoiseAdj')
               KSpace{iMeas}.(typeStr).data = [];
            end
            if Control.User.readRTFeedback & strcmp(typeStr,'RTFeedback')
               KSpace{iMeas}.(typeStr).data = [];
            end
         end
         progressPercent = progressPercent + 0.25;
         if ~(MrParcRaidFileCell{iMeas}.Control.User.noGui || MrParcRaidFileCell{iMeas}.Control.User.silent)
            waitbar(progressPercent, Control.GUI.hWaitBar, 'Allocating arrays');
         end
      end

      iScanRead = uint64(0);
      if ~(MrParcRaidFileCell{iMeas}.Control.User.noGui || MrParcRaidFileCell{iMeas}.Control.User.silent)
         time0 = tic;
         time00 = tic;
      end

      for iBlock=1:length(MdhBlock)

         % Skip SYNCDATA
         if MdhBlock(iBlock).isSyncData || MdhBlock(iBlock).isAcqEnd
            continue
         end
         iReadout = uint64(1):uint64(MdhBlock(iBlock).nSamplesInScan);

         maxScanPerRead = max(idivide(MAX_BYTES_PER_FREAD,MdhBlock(iBlock).length,'fix'),uint64(1));
         nScanToRead = min(maxScanPerRead,MdhBlock(iBlock).nScan);

         fseek(fid,MdhBlock(iBlock).offsetInFile,'bof');
         nScanRemainInBlock = MdhBlock(iBlock).nScan;

         while nScanRemainInBlock > 0

            if ~(MrParcRaidFileCell{iMeas}.Control.User.noGui || MrParcRaidFileCell{iMeas}.Control.User.silent)
               dTimeElapsed = toc(time0);
               rate = dTimeElapsed / (double(iScanRead) + 1);
               timeRemain = double(nScanInMeas - iScanRead) * rate;
               textMessage = sprintf('Pass 3 of 3\nReading scans\nTime remaining (s): %d', uint32(timeRemain));
               waitbar(double(iScanRead)/double(nScanInMeas), Control.GUI.hWaitBar, textMessage);
            end
            blockOfScans = fread(fid,[MdhBlock(iBlock).length nScanToRead],'uchar=>uchar');
            MdhParam = extract_mdh_parameters(blockOfScans, 'scan', {'timeStamp', 'cutOffDataPre', 'cutOffDataPost', 'indexLIN','indexACQ','indexSLC','indexPAR','indexECO','indexPHS','indexREP','indexSET','indexSEG'});
            MdhBitFlags = extract_mdh_bitflags(blockOfScans, {'isAcqEnd', 'isNoiseAdj', 'isPhaseCor', 'isRTFeedback', 'isRawDataCorrection','isReflect'});

            for indexStr = {'indexLIN','indexACQ','indexSLC','indexPAR','indexECO','indexPHS','indexREP','indexSET','indexSEG'}
               MdhParam.(indexStr{1}) = uint64(MdhParam.(indexStr{1}));
            end
            channelID = uint16(0) * zeros(MdhBlock(iBlock).nChannelUsed,nScanToRead,'uint16');
            for iCha = 1:MdhBlock(iBlock).nChannelUsed
               offset  = uint64(MDH_SCANSIZE) + uint64(iCha-1) * sizeChanInBytes;
               channelID(iCha,:) = reshape(extract_single_mdh_parameter(blockOfScans(offset+1:offset+sizeChanInBytes,:),'chan','channelID'),[1,nScanToRead]);
            end
            sizeSamplesInBytes = MdhBlock(iBlock).nSamplesInScan * uint64(8);
            for iCha = 1:MdhBlock(iBlock).nChannelUsed

               channelIDOut = squeeze(channelID(iCha,:));
               offset  = uint64(MDH_SCANSIZE) + uint64(iCha-1) * sizeChanInBytes + MDH_CHANSIZE + uint64(1);

               tempArr = blockOfScans(offset:offset+sizeSamplesInBytes-1,:);
               tempArr = reshape(typecast(tempArr(:),'single'), [2*MdhBlock(iBlock).nSamplesInScan nScanToRead]);
               tempArr = complex(tempArr(1:2:end,:),tempArr(2:2:end,:));

               % Handle all the reflections at once
               indexReflect = find(MdhBitFlags.isReflect);
               if ~isempty(indexReflect)
                  tempArr(:,indexReflect) = flipud(tempArr(:,indexReflect));
               end

               % Handle Pre-cutoff
               indexPre = find(MdhParam.cutOffDataPre > 0 & ~MdhBitFlags.isNoiseAdj);
               if ~isempty(indexPre)
                  cutOffDataPreUniq = unique(MdhParam.cutOffDataPre(indexPre));
                  for iUniq = 1:numel(cutOffDataPreUniq)
                     indexCut = find(MdhParam.cutOffDataPre(indexPre) == cutOffDataPreUniq(iUniq));
                     tempArr(1:cutOffDataPreUniq(iUniq),indexPre(indexCut)) = ...
                        complex(zeros(cutOffDataPreUniq(iUniq),numel(indexCut),'single'));
                  end % loop over uniq cuts
               end

               % Handle Post-cutoff
               indexPost = find(MdhParam.cutOffDataPost > 0 & ~MdhBitFlags.isNoiseAdj);
               if ~isempty(indexPost)
                  cutOffDataPostUniq = unique(MdhParam.cutOffDataPost(indexPost));
                  for iUniq = 1:numel(cutOffDataPostUniq)
                     indexCut = find(MdhParam.cutOffDataPost(indexPost) == cutOffDataPostUniq(iUniq));
                     tempArr(1:cutOffDataPostUniq(iUniq),indexPost(indexCut)) = ...
                        complex(zeros(cutOffDataPostUniq(iUniq),numel(indexCut),'single'));
                  end % loop over uniq cuts
               end

               maskNormal = ~MdhBitFlags.isAcqEnd;
               for typeCell = {'NoiseAdj' 'PhaseCor' 'RTFeedback' 'Normal'}
                  typeStr = typeCell{1};
                  bitField = strcat('is',typeStr);
                  readField = strcat('read',typeStr);
                  flagNormal = strcmp(typeStr,'Normal');
                  flagPhaseCor = strcmp(typeStr,'PhaseCor');
                  flagNSEG = n.SEG.(typeStr) > 1;
                  if ~flagNormal
                     maskNormal = maskNormal & (MdhBitFlags.(bitField) == 0);
                     index.(typeStr) = find(MdhBitFlags.(bitField) & ~MdhBitFlags.isAcqEnd);
                  else
                     index.(typeStr) = find(maskNormal);
                  end
                  if ~(MrParcRaidFileCell{iMeas}.Control.User.(readField) && (numel(index.(typeStr)) > 0))
                     continue;
                  end
                  % Pass through to sort optimally
                  kScanTest = zeros(numel(index.(typeStr)),1,'uint64');
                  for kScan = 1:numel(index.(typeStr))
                     jScan = index.(typeStr)(kScan);
                     if flagNormal
                        keyCHACoilSelectMeas = channelIDOut(jScan) + 1;
                        keyCHACoilSelect = CoilSelectMeasMap(keyCHACoilSelectMeas).tElement;
                        indexCHA = find(channelIDUniq == (keyCHACoilSelectMeas-1)) - 1;
                     else
                        indexCHA = find(channelIDUniq == channelIDOut(jScan)) - 1;
                     end
                     indexOutTemp = uint64([indexCHA MdhParam.indexLIN(jScan) MdhParam.indexSLC(jScan) MdhParam.indexPAR(jScan) MdhParam.indexACQ(jScan) MdhParam.indexECO(jScan) MdhParam.indexPHS(jScan) MdhParam.indexREP(jScan) MdhParam.indexSET(jScan)]);
                     indexOut.(typeStr) = mod(indexOutTemp(order.(typeStr)),dimOut.(typeStr));
                     kScanTest(kScan,1) = sum(indexOut.(typeStr)(:) .* dimOutProd.(typeStr)(:)) + uint64(1);
                  end
                  [kScanSorted,kScanSortIndex] = sort(kScanTest);
                  for kScan = 1:numel(index.(typeStr))
                     if ~(MrParcRaidFileCell{iMeas}.Control.User.noGui || MrParcRaidFileCell{iMeas}.Control.User.silent)
                        dTimeElapsed00 = toc(time00);
                        if dTimeElapsed00 > 5
                           dTimeElapsed = toc(time0);
                           rate = dTimeElapsed / (double(iScanRead+kScan) + 1);
                           timeRemain = double(nScanInMeas - iScanRead - kScan) * rate;
                           time00 = tic;
                           textMessage = sprintf('Pass 3 of 3\nReading scans\nTime remaining (s): %d', uint32(timeRemain));
                           waitbar(double(iScanRead)/double(nScanInMeas), Control.GUI.hWaitBar, textMessage);
                        end
                        rate = dTimeElapsed / (double(iScanRead) + 1);
                        timeRemain = double(nScanInMeas - iScanRead) * rate;
                     end
                     jScan = index.(typeStr)(kScanSortIndex(kScan));
                     if flagNormal
                        keyCHACoilSelectMeas = channelIDOut(jScan) + 1;
                        keyCHACoilSelect = CoilSelectMeasMap(keyCHACoilSelectMeas).tElement;
                        indexCHA = find(channelIDUniq == (keyCHACoilSelectMeas-1)) - 1;
                     else
                        indexCHA = find(channelIDUniq == channelIDOut(jScan)) - 1;
                     end
                     indexOutTemp = uint64([indexCHA MdhParam.indexLIN(jScan) MdhParam.indexSLC(jScan) MdhParam.indexPAR(jScan) MdhParam.indexACQ(jScan) MdhParam.indexECO(jScan) MdhParam.indexPHS(jScan) MdhParam.indexREP(jScan) MdhParam.indexSET(jScan)]);
                     temp = tempArr(:,jScan);
                     indexOut.(typeStr) = mod(indexOutTemp(order.(typeStr)),dimOut.(typeStr));
                     iScan.(typeStr) = sum(indexOut.(typeStr)(:) .* dimOutProd.(typeStr)(:)) + uint64(1);
                     if flagPhaseCor
                        KSpace{iMeas}.(typeStr).data(iReadout,iScan.(typeStr),MdhParam.indexSEG(jScan)+1) = temp;
                        if Control.User.readTimeStamp
                           KSpace{iMeas}.(typeStr).TimeStamp.count(iScan.(typeStr),MdhParam.indexSEG(jScan)+1) = KSpace{iMeas}.(typeStr).TimeStamp.count(iScan.(typeStr),MdhParam.indexSEG(jScan)+1) +1;
                           nTimeStamp = KSpace{iMeas}.(typeStr).TimeStamp.count(iScan.(typeStr),MdhParam.indexSEG(jScan)+1); 
                           KSpace{iMeas}.(typeStr).TimeStamp.value(iScan.(typeStr),MdhParam.indexSEG(jScan)+1,nTimeStamp) = MdhParam.timeStamp(jScan);
                        end
                     else
                        if ~flagNormal
                           KSpace{iMeas}.(typeStr).data(iReadout,iScan.(typeStr)) = temp;
                        else
                           if MdhBitFlags.isRawDataCorrection(jScan) && Control.User.applyRawDataCorrection
                              temp = CoilSelectMap(keyCHACoilSelect).rawDataCorrectionFactor*temp;
                           end
                           KSpace{iMeas}.Normal.data(iReadout,iScan.Normal) = CoilSelectMeasMap(keyCHACoilSelectMeas).flFFTCorrectionFactor*temp;
                        end
                        if flagNSEG
                           KSpace{iMeas}.(typeStr).SEG(1,iScan.(typeStr)) = MdhParam.indexSEG(jScan) + 1;
                        end
                        if Control.User.readTimeStamp
                           KSpace{iMeas}.(typeStr).TimeStamp.count(iScan.(typeStr)) = KSpace{iMeas}.(typeStr).TimeStamp.count(iScan.(typeStr)) +1;
                           nTimeStamp = KSpace{iMeas}.(typeStr).TimeStamp.count(iScan.(typeStr)); 
                           KSpace{iMeas}.(typeStr).TimeStamp.value(iScan.(typeStr),nTimeStamp) = MdhParam.timeStamp(jScan);
                        end
                     end
                  end
               end

            end % endfor over channels

            iScanRead = iScanRead + nScanToRead;
            nScanRemainInBlock = MdhBlock(iBlock).nScan - iScanRead;
            nScanToRead = min(nScanToRead,nScanRemainInBlock);

         end % while over sub-blocks of lines
      end % for over MdhBlocks
      clear blockOfScans;
      clear index;
      clear tempArr;
      clear channelID channelIDOut MdhParam MdhBitFlags;
      clear iScanRead nScanToRead;
      if ~(MrParcRaidFileCell{iMeas}.Control.User.noGui || MrParcRaidFileCell{iMeas}.Control.User.silent)
         waitbar(0.5, Control.GUI.hWaitBar, 'All lines read from file.  Reshaping.');
      end

      % Reshape the results
      for typeCell = {'NoiseAdj' 'PhaseCor' 'RTFeedback' 'Normal'}
         typeStr = typeCell{1};
         readField = strcat('read',typeStr);
         if ~MrParcRaidFileCell{iMeas}.Control.User.(readField) 
            continue;
         end
         flagPhaseCor = strcmp(typeStr,'PhaseCor');
         flagNSEG = n.SEG.(typeStr) > 1;
         indexDimNonZero.(typeStr) = find(dimOut.(typeStr) > 1);
         if isempty(indexDimNonZero.(typeStr))
            indexDimNonZeroTemp = find(dimOut.(typeStr)(1:end) > 0); 
            indexDimNonZero.(typeStr) = zeros(1,1);
            indexDimNonZero.(typeStr)(1) = indexDimNonZeroTemp(1);
         end
         if flagPhaseCor
            if ~isempty(KSpace{iMeas}.(typeStr).data)
               KSpace{iMeas}.PhaseCor.data = reshape(KSpace{iMeas}.PhaseCor.data, [n.COL dimOut.PhaseCor(indexDimNonZero.PhaseCor) n.SEG.PhaseCor]);
               if Control.User.readTimeStamp
                  dSize = max(KSpace{iMeas}.PhaseCor.TimeStamp.count(:));
                  KSpace{iMeas}.PhaseCor.TimeStamp.count = reshape(KSpace{iMeas}.PhaseCor.TimeStamp.count, [dimOut.PhaseCor(indexDimNonZero.PhaseCor) n.SEG.PhaseCor]);
                  KSpace{iMeas}.PhaseCor.TimeStamp.value = reshape(KSpace{iMeas}.PhaseCor.TimeStamp.value, [dimOut.PhaseCor(indexDimNonZero.PhaseCor) n.SEG.PhaseCor dSize]);
               end
            end
         else
            if ~isempty(KSpace{iMeas}.(typeStr).data)
               KSpace{iMeas}.(typeStr).data = reshape(KSpace{iMeas}.(typeStr).data, [n.COL dimOut.(typeStr)(indexDimNonZero.(typeStr))]);
               if flagNSEG
                  KSpace{iMeas}.(typeStr).SEG = squeeze(reshape(KSpace{iMeas}.(typeStr).SEG, [1 dimOut.(typeStr)(indexDimNonZero.(typeStr))]));
               end
               if Control.User.readTimeStamp
                  dSize = max(KSpace{iMeas}.(typeStr).TimeStamp.count(:));
                  KSpace{iMeas}.(typeStr).TimeStamp.count = reshape(KSpace{iMeas}.(typeStr).TimeStamp.count, dimOut.(typeStr)(indexDimNonZero.(typeStr)));
                  KSpace{iMeas}.(typeStr).TimeStamp.value = reshape(KSpace{iMeas}.(typeStr).TimeStamp.value, [dimOut.(typeStr)(indexDimNonZero.(typeStr)) dSize]);
               end
            end
         end
      end

      if MrParcRaidFileCell{iMeas}.Control.User.shiftDCToMatrixCenter == true
         for typeCell = {'NoiseAdj' 'PhaseCor' 'RTFeedback' 'Normal'}
            typeStr = typeCell{1};
            flagPhaseCor = strcmp(typeStr,'PhaseCor');
            flagNSEG = n.SEG.(typeStr) > 1;
            readField = strcat('read',typeStr);
            if MrParcRaidFileCell{iMeas}.Control.User.(readField) && (nScan.(typeStr) > 0)
               allShift = dimOut.(typeStr) * 0;
               colShift = idivide(n.COL - nSamplesInScan, uint64(2), 'fix');
               if dimOut.(typeStr)(2) > 1
                  linShift = uint64(dimOut.(typeStr)(2)/2) - uint64(kSpaceCentreLineNo);
                  allShift(2) = linShift;
               end
               if dimOut.(typeStr)(4) > 1
                  parShift = uint64(dimOut.(typeStr)(4)/2) - uint64(kSpaceCentrePartitionNo);
                  allShift(4) = parShift;
               end
               allShift = allShift(order.(typeStr));
               if flagPhaseCor
                  allShift = [colShift allShift(indexDimNonZero.(typeStr)) 0];
               else
                  allShift = [colShift allShift(indexDimNonZero.(typeStr))];
               end
               KSpace{iMeas}.(typeStr).data = circshift(KSpace{iMeas}.(typeStr).data,allShift);
            end
         end
      end

      if ~(MrParcRaidFileCell{iMeas}.Control.User.noGui || MrParcRaidFileCell{iMeas}.Control.User.silent)
         waitbar(1.0, Control.GUI.hWaitBar, 'KSpace complete.');
      end

      tMsg = ['No normal scans read this measUID: ' MrParcRaidFileCell{iMeas}.Control.File.measUIDString];
      if nScan.Normal > 0
         tMsg = sprintf('Normal KSpace dimensions:\n   COL=%d', n.COL);
         dimNamesOut.Normal = dimNames(order.Normal);
         for iDim=1:numel(indexDimNonZero.Normal)
            tMsg = strcat(tMsg, sprintf('\n   %s=%d', ...
                                        dimNamesOut.Normal{indexDimNonZero.Normal(iDim)}, ...
                                        dimOut.Normal(indexDimNonZero.Normal(iDim))));
         end
         if (n.SEG.Normal > 1) 
            tMsg = strcat(tMsg, sprintf('\n   SEG=%d', n.SEG.Normal));
         end
         KSpace{iMeas}.Normal.dim = ['COL', dimNamesOut.Normal(indexDimNonZero.Normal)];
      end

      if nScan.NoiseAdj > 0
         tMsg = sprintf('NoiseAdj KSpace dimensions:\n   COL=%d', n.COL);
         dimNamesOut.NoiseAdj = dimNames(order.NoiseAdj);
         for iDim=1:numel(indexDimNonZero.NoiseAdj)
            tMsg = strcat(tMsg, sprintf('\n   %s=%d', ...
                                        dimNamesOut.NoiseAdj{indexDimNonZero.NoiseAdj(iDim)}, ...
                                        dimOut.NoiseAdj(indexDimNonZero.NoiseAdj(iDim))));
         end
         if (n.SEG.NoiseAdj > 1) 
            tMsg = strcat(tMsg, sprintf('\n   SEG=%d', n.SEG.NoiseAdj));
         end
         KSpace{iMeas}.NoiseAdj.dim = ['COL', dimNamesOut.NoiseAdj(indexDimNonZero.NoiseAdj)];
      end

      if nScan.RTFeedback > 0
         tMsg = sprintf('RTFeedback KSpace dimensions:\n   COL=%d', n.COL);
         dimNamesOut.RTFeedback = dimNames(order.RTFeedback);
         for iDim=1:numel(indexDimNonZero.RTFeedback)
            tMsg = strcat(tMsg, sprintf('\n   %s=%d', ...
                                        dimNamesOut.RTFeedback{indexDimNonZero.RTFeedback(iDim)}, ...
                                        dimOut.RTFeedback(indexDimNonZero.RTFeedback(iDim))));
         end
         if (n.SEG.RTFeedback > 1) 
            tMsg = strcat(tMsg, sprintf('\n   SEG=%d', n.SEG.RTFeedback));
         end
         KSpace{iMeas}.RTFeedback.dim = ['COL', dimNamesOut.RTFeedback(indexDimNonZero.RTFeedback)];
      end

      if nScan.PhaseCor > 0
         tMsg = strcat(tMsg, sprintf('\nPhase Cor dimensions:\n   col=%d', n.COL));
         dimNamesOut.PhaseCor = dimNames(order.PhaseCor);
         for iDim=1:numel(indexDimNonZero.PhaseCor)
            tMsg = strcat(tMsg, sprintf('\n   %s=%d', ...
                                        dimNamesOut.PhaseCor{indexDimNonZero.PhaseCor(iDim)}, ...
                                        dimOut.PhaseCor(indexDimNonZero.PhaseCor(iDim))));
         end
         tMsg = strcat(tMsg, sprintf('\n   SEG=%d', n.SEG.PhaseCor));
         KSpace{iMeas}.PhaseCor.dim = ['COL', dimNamesOut.PhaseCor(indexDimNonZero.PhaseCor), 'SEG'];
      end

      if ~MrParcRaidFileCell{iMeas}.Control.User.silent
         display_message(MrParcRaidFileCell{iMeas}.Control.User.noGui, ...
                         'WAIT', tMsg,'Output KSpace Dimensions');
      end

      KSpace{iMeas}.Normal.Protocol.Dim = Dim;
      KSpace{iMeas}.Normal.Protocol.CoilSelectMap = CoilSelectMap;
      KSpace{iMeas}.Normal.Protocol.CoilSelectMeasMap = CoilSelectMeasMap;

      clear Dim CoilSelectMap CoilSelectMeasMap n;

   end % endfor over all measurements in file
   if ~(MrParcRaidFileCell{iMeas}.Control.User.noGui || MrParcRaidFileCell{iMeas}.Control.User.silent)
      close(Control.GUI.hWaitBar);
   end
   fclose(fid);

   if Control.User.readPMU
      varargout{Control.User.iArgOutPMU} = PMUOutCell;
   end
   if and(Control.User.foundWiPMemBlock,Control.User.readWiPMemBlock)
      varargout{Control.User.iArgOutWiPMemBlock} = WiPMemBlock;
   end
   if and(Control.User.foundRelSliceNumber,Control.User.readRelSliceNumber)
      varargout{Control.User.iArgOutRelSliceNumber} = RelativeSliceNumber;
   end

end

% =============================================================================
function UserInputDefault = define_user_input_default_struct_array
% -----------------------------------------------------------------------------
% Return struct array of user defaults:
%    UserInputDefault(i).varArgIn - string to expect from varargin
%    UserInputDefault(i).fieldName - defines structure field name to associate 
%                                with input.  If empty string, then same as 
%                                varArgIn with first letter lower case.
%    UserInputDefault(i).defaultValue - default value
%    UserInputDefault(i).type - currently: 'logical','uint64'
%    UserInputDefault(i).isArgOut - boolean
%    UserInputDefault(i).description
% -----------------------------------------------------------------------------

   userInputFieldNames = ...
         {{'varArgIn'},               {'fieldName'}, {'defaultValue'}, {'type'},    {'isArgOut'}, {'description'} };
   userInputValues = ...
      { ...
         {{'ApplyRawDataCorrection'}, {''},               {false},     {'logical'}, {false},      {'Apply raw data correction to (usually central TSE) lines based on mdh flag and gain parameters in the header.'}}, ...
         {{'LastMeasDatOnly'},        {''},               {false},     {'logical'}, {false},      {'Read only the last meas.dat block in the multi-measure meas.dat file.'}}, ...
         {{'NoFillToFullFourier'},    {''},               {false},     {'logical'}, {false},      {'Do NOT Fill k-space matrix out to complete FFT size in case of partial-Fourier acquisitions.'}}, ...
         {{'NoGui'},                  {''},               {false},     {'logical'}, {false},      {'Do not display GUI elements.'}}, ...
         {{'ReadNoiseAdj'},           {''},               {false},     {'logical'}, {false},      {'Read noise adjust lines.  Returns kspace.  User must supply output variable.'}}, ...
         {{'ReadPhaseCor'},           {''},               {false},     {'logical'}, {false},      {'Read phase correction lines.  Returns kspace.  User must supply output variable.'}}, ...
         {{'ReadPMU'},                {''},               {false},     {'logical'}, {true},       {'Reads in the embedded PMU times and waveforms from the meas.dat header.  Usually requires gated sequences.'}} ...
         {{'ReadRelSliceNumber'},     {''},               {false},     {'logical'}, {true},       {'Read relative slice number from header.  Necessary to sort 2D acquisitions.'}} ...
         {{'ReadRTFeedback'},         {''},               {false},     {'logical'}, {false},      {'Read RT feedback lines.  Returns kspace.  User must supply output variable.'}}, ...
         {{'ReadTimeStamp'},          {''},               {false},     {'logical'}, {false},      {'Adds additional field to KSpace structures storing time stamps corresponding to KSpace lines.'}}, ...
         {{'ReadWiPMemBlock'},        {''},               {false},     {'logical'}, {true},       {'Extract WiPMemBlock from header.  User must provide return variable in var arg out.'}}, ...
         {{'ResetFFTScale'},          {''},               {false},     {'logical'}, {false},      {'Reset scale for each coil (read from FFT Correction Factors) to one.'}}, ...
         {{'ShiftDCToMatrixCenter'},  {''},               {false},     {'logical'}, {false},      {'Shifts output k-space such that DC occurs at n/2+1 for all dimensions.'}}, ...
         {{'Silent'},                 {''},               {false},     {'logical'}, {false},      {'Turns off all GUI and output messages.  User must supply all necessary arguments.'}}, ...
         {{'SwapSegAndEco'},          {''},               {false},     {'logical'}, {false},      {'Swap incoming segment and echo indices.  Hack to permit reading of HIFU navigators encoded with SEG indices.'}}, ...
      };

   for iUserInput = 1:length(userInputValues)
      for iField = 1:length(userInputFieldNames)
         UserInputDefault(iUserInput).(userInputFieldNames{iField}{1}) = userInputValues{iUserInput}{iField}{1};
      end     
   end

end

% =============================================================================
function UserInputDefault = derive_user_input_default_fieldname(UserInputDefault)
% -----------------------------------------------------------------------------
% Derives fieldnames to used later for User Input struct.  
% Names based on variation of possible user inputs.
%    ex.  User input 'ReadSomething'  becomes field name 'readSomething'
% -----------------------------------------------------------------------------

   % Setup fieldNames
   for iUserInput = 1:length(UserInputDefault)
      if isempty(UserInputDefault(iUserInput).fieldName)
         fieldName = UserInputDefault(iUserInput).varArgIn;
         fieldName(1) = lower(fieldName(1));
         UserInputDefault(iUserInput).fieldName = fieldName;
      end
   end

end

% =============================================================================
function User = convert_user_input_default_to_struct(UserInputDefault)
% -----------------------------------------------------------------------------
% Set up user input structure with default values.
% -----------------------------------------------------------------------------

   for iUserInput = 1:length(UserInputDefault)
      User.(UserInputDefault(iUserInput).fieldName) = UserInputDefault(iUserInput).defaultValue;
      fieldName = ['found',UserInputDefault(iUserInput).varArgIn(5:end)];
      User.(fieldName) = false;
   end

end

% =============================================================================
function filenameIn = check_filename_input(filenameIn, nArgIn, Control)
% -----------------------------------------------------------------------------
% Check filenameIn and query if empty.
% -----------------------------------------------------------------------------

   if nArgIn == 0
      if ~Control.User.noGui
         filenameIn = uigetfile('*.dat','Select File to Read');
      else
         if ~Control.User.silent
            textMessage = 'Must provide file name to read';
            display_message(Control.User.noGui,'ERROR',textMessage);
         end
         return
      end
   end

end

% =============================================================================
function filenameInExists = check_filename_exists(filenameIn, Control)
% -----------------------------------------------------------------------------
% Check for file existence.
% -----------------------------------------------------------------------------
   filenameInExists = (exist(filenameIn, 'file') == 2);
   if ~filenameInExists
      if ~Control.User.silent
         textMessage = ['File does not exist: ' filenameIn];
         display_message(Control.User.noGui,'ERROR',textMessage,'Error');
      end
   end
end

% =============================================================================
function Control = modify_with_inputs_outputs(optionsIn, Control, UserInputDefault)
% -----------------------------------------------------------------------------
% Adjust controls based on user inputs/outputs
% -----------------------------------------------------------------------------

   Control.nArgOut = 0;

   iOption = 1;
   while iOption <= length(optionsIn)

      flagFoundMatch = false;
      for iUserInput = 1:length(UserInputDefault)
         flagFoundMatch = flagFoundMatch || strcmpi(optionsIn{iOption},UserInputDefault(iUserInput).varArgIn);
         if flagFoundMatch
            break;
         end
      end
      if flagFoundMatch
         fieldName = UserInputDefault(iUserInput).fieldName;
         if strcmp(UserInputDefault(iUserInput).type, 'logical')
            Control.User.(fieldName) = true;
            if UserInputDefault(iUserInput).isArgOut
               fieldName = ['iArgOut',UserInputDefault(iUserInput).varArgIn(5:end)];
               Control.nArgOut = Control.nArgOut + 1;
               Control.User.(fieldName) = Control.nArgOut;
            end
         end
         if strcmp(UserInputDefault(iUserInput).type, 'uint64')
            iOption = iOption + 1;
            Control.User.(fieldName) = optionsIn{iOption}(1);
         end
      end
      iOption = iOption + 1;
   end % while over options

end

% =============================================================================
function goodToGo = are_inputs_and_outputs_consistent(nOptionsIn, nArgOutMain, Control)
% -----------------------------------------------------------------------------
% Check for consistentcy between argin and argout when a given argin might
%    lead to more data in the output.
% -----------------------------------------------------------------------------

   goodToGo = true;
   if (nOptionsIn > 0) && ((nArgOutMain-1) ~= Control.nArgOut)
      goodToGo = false;
      if ~Control.User.silent
         textMessage = 'You must have an equal number of kSpace output variables for all the kSpace you expect to read.'; 
         display_message(Control.User.noGui,'ERROR',textMessage,'Insufficient Outputs');
      end
   end

   if Control.User.shiftDCToMatrixCenter && Control.User.readTimeStamp
      goodToGo = false;
      if ~Control.User.silent
         textMessage = 'You cannot turn on ShiftDCToMatrixCenter when using ReadTimeStamp';
         display_message(Control.User.noGui,'ERROR',textMessage,'Insufficient Outputs');
      end
   end

end

% =============================================================================
function Control = get_controls_from_user_input(filenameIn, optionsIn, numArgIn, numArgOut)
% -----------------------------------------------------------------------------
% Setup main control structure using user inputs.
% -----------------------------------------------------------------------------

   Control.GUI.hWaitBar = []; % Initialize to blank

   nOptionsIn = numArgIn - 1;

   UserInputDefault = define_user_input_default_struct_array;
   UserInputDefault = derive_user_input_default_fieldname(UserInputDefault);

   Control.User = convert_user_input_default_to_struct(UserInputDefault);
   Control.User.readNormal = 1; % Always
   Control = modify_with_inputs_outputs(optionsIn, Control, UserInputDefault);

   Control.File.In.name = check_filename_input(filenameIn, nargin, Control);
   if ~check_filename_exists(Control.File.In.name, Control)
      Control = [];
      return
   end

   if ~are_inputs_and_outputs_consistent(nOptionsIn, numArgOut, Control)
      Control = [];
      return 
   end

   if ~Control.User.silent
      print_setup_status(Control.File.In.name, UserInputDefault, Control);
   end

end

% =============================================================================
function [isVD13, MrParcRaidFileHeader, MrParcRaidFileCell] = is_vd13_file(Control, fid)
% -----------------------------------------------------------------------------
% Check that this is a VD13 file.
% -----------------------------------------------------------------------------

   MrParcRaidFileHeader = [];
   MrParcRaidFileCell = {};
   isVD13 = 0;
   FileHeader.totalSize  = fread(fid,1,'uint32=>uint32');
   FileHeader.nMeas = fread(fid,1,'uint32=>uint32');
   if not(and(FileHeader.totalSize < 10000, FileHeader.nMeas <= 64))
      if ~Control.User.silent
         display_message(Control.User.noGui,'ERROR',['Appears to be VB data: ' Control.File.In.name],'Error');
      end
      return 
   else
      MrParcRaidFileHeader = FileHeader;
      isVD13 = 1;
   end

   % Extract measurement parameters
   for iMeas = 1:MrParcRaidFileHeader.nMeas

      MrParcRaidFileCell{iMeas}.measID = fread(fid,1,'uint32=>uint32');
      MrParcRaidFileCell{iMeas}.fileID = fread(fid,1,'uint32=>uint32');
      MrParcRaidFileCell{iMeas}.offsetInFile = fread(fid,1,'uint64=>uint64');
      MrParcRaidFileCell{iMeas}.totalSize = fread(fid,1,'uint64=>uint64');
      MrParcRaidFileCell{iMeas}.patName = fread(fid,64,'uchar=>char');
      MrParcRaidFileCell{iMeas}.protName = fread(fid,64,'uchar=>char');
      MrParcRaidFileCell{iMeas}.Control = Control;

   end

end

% =============================================================================
function [Control, MrParcRaidFileCell, Dim, CoilSelectMap, ...
          CoilSelectMeasMap, RelativeSliceNumber, WiPMemBlock] = ...
             parse_meas_text_header(Control,fid,MrParcRaidFileCell,iMeas);
% -----------------------------------------------------------------------------
% Parse the text portion of the measurement to extract required dimensions,
%    scaling factors, slice ordering, etc.
% -----------------------------------------------------------------------------

   % Skip to start of this measurement
   fseek(fid,MrParcRaidFileCell{iMeas}.offsetInFile,'bof');

   % TextHeader extraction
   Dim.dMeasHeaderSize = fread(fid,1,'int32');
   Dim.ullMeasHeaderSize = uint64(Dim.dMeasHeaderSize);
   textHeader = fread(fid,Dim.dMeasHeaderSize-4,'uchar=>char');
   textHeader = textHeader';
   TextHeader.complete = textHeader;
   clear textHeader;

   findThisStart = 'MeasYaps';
   findThisEnd = 'Phoenix';
   pStart = strfind(TextHeader.complete,findThisStart) + length(findThisStart) + 5;
   pEnd = strfind(TextHeader.complete,findThisEnd) - 3;
   TextHeader.measYapsAscConv = TextHeader.complete(pStart:pEnd);

   searchCells = {{'Dim.Raw.readoutOversampleFactor'  'complete' '.*flReadoutOSFactor[^\s]*\s*{\s*<Precision>\s*[0-9]*\s*([0-9.]+)' 'str2num' 2.0} ...
                  {'Dim.Raw.nCha' 'measYapsAscConv' '\[0\]\.asList\[([0-9]+)\]\.lRxChannelConnected\s*=\s*'       'length'              } ...
                  {'Dim.Recon.nFourierPartitions'  'complete'  '.*iNoOfFourierPartitions[^\s]*\s*{\s*([0-9.]+)'   'str2num'            1} ...
                  {'Dim.Recon.phaseEncodeFTLength'  'complete'      '.*iPEFTLength[^\s]*\s*{\s*([0-9.]+)'         'str2num'             } ...
                  {'Dim.Recon.partitionFTLength'  'complete'        '.*i3DFTLength[^\s]*\s*{\s*([0-9.]+)'         'str2num'             } ...
                  {'Control.File.scanDimension'   'measYapsAscConv' 'sKSpace.ucDimension[^=]+=\s*([0-9]+)\s*'     'str2num'             } ...
                  {'Dim.Recon.nLin'               'complete'        'iPEFTLen[^{}]+{\s*([0-9]+)\s*}'              'str2double'          } ...
                  {'Control.File.measUIDString'   'complete'        'MeasUID[^{}]+{\s*([0-9]+)\s*}'               'string'              } ...
                  {'Dim.Recon.nPar'               'measYapsAscConv' 'sKSpace.lPartitions[^=]+=\s*([0-9]+)\s*'     'str2num'             } ...
                  {'Dim.Recon.nSlc'               'measYapsAscConv' 'sSliceArray.lSize[^=]+=\s*([0-9]+)\s*'       'str2num'             } ...
                  {'Dim.Raw.nEco'                 'measYapsAscConv' 'lContrasts\s*=\s*([0-9]+)\s*'                'str2num'            1} ...
                  {'Dim.Raw.nSet'                 'measYapsAscConv' 'lSets\s*=\s*([0-9]+)\s*'                     'str2num'            1} ...
                  {'Dim.Raw.nAcq'                 'measYapsAscConv' 'lAverages\s*=\s*([0-9]+)\s*'                 'str2num'            1} ...
                  {'Dim.Raw.nRep'                 'measYapsAscConv' 'lRepetitions\s*=\s*([0-9]+)\s*'              'str2num'            1} ...
                  {'Dim.Raw.nPhs'                 'measYapsAscConv' 'sPhysioImaging.lPhases\s*=\s*([0-9]+)\s*'    'str2num'            1} ...
                  {'Control.File.epiFactor'       'measYapsAscConv' 'sFastImaging.lEPIFactor\s*=\s*([0-9]+)\s*'   'str2num'             } ...
                  {'Control.File.turboFactor'     'measYapsAscConv' 'sFastImaging.lTurboFactor\s*=\s*([0-9]+)\s*' 'str2num'            1} ...
                  {'Dim.Raw.nMaxRxChannels'       'complete'        '.*iMaxNoOfRxChannels[^\s]*\s*{\s*([0-9]+)'   'str2num'          128} ...
                  {'Dim.Raw.nPhaseCor'            'complete'        '.*lNoOfPhaseCorrScans[^\s]*\s*{\s*([0-9]+)'  'str2num'            0} ...
                  {'Dim.Recon.nFourierColumns'    'complete'        '.*iNoOfFourierColumns[^\s]*\s*{\s*([0-9]+)'  'str2num'             } ...
                  {'Dim.Recon.nFourierLines'      'complete'        '.*NoOfFourierLines[^\s]*\s*{\s*([0-9.]+)'    'str2num'             } ...
                  {'Dim.Recon.readoutFTLength'    'complete'        '.*iRoFTLength[^\s]*\s*{\s*([0-9.]+)'         'str2num'             }};
   for searchCell = searchCells
      flagHasDefault = length(searchCell{1}) > 4;
      address = strsplit(searchCell{1}{1},'.');
      searchField = searchCell{1}{2};
      searchString = searchCell{1}{3};
      conversion = searchCell{1}{4};
      if flagHasDefault
         defaultValue = searchCell{1}{5};
      end
      searchResult= regexp(TextHeader.(searchField),searchString,'tokens');
      if ~isempty(searchResult)
         if strcmp(conversion,'str2num')
            searchValue = str2num(searchResult{1}{1});
         elseif strcmp(conversion,'string')
            searchValue = searchResult{1}{1};
         elseif strcmp(conversion,'length')
            searchValue = length(searchResult);
         elseif strcmp(conversion,'str2double')
            searchValue = str2double(searchResult{1}{1});
         else
            Control = [];
            MrParcRaidFileCell = [];
            Dim = [];
            CoilSelectMap = [];
            CoilSelectMeasMap = [];
            RelativeSliceNumber = [];
            WiPMemBlock = [];
          return
         end
      elseif flagHasDefault
         searchValue = defaultValue;
      else
         continue;
      end
      if strcmp(address{1},'Dim')
         Dim.(address{2}).(address{3}) = searchValue;
      end
      if strcmp(address{1},'Control')
         MrParcRaidFileCell{iMeas}.Control.(address{2}).(address{3}) = searchValue;
      end
   end

   if MrParcRaidFileCell{iMeas}.Control.File.scanDimension == 4
      MrParcRaidFileCell{iMeas}.Control.File.is3D = true;
   else
      MrParcRaidFileCell{iMeas}.Control.File.is3D = false;
   end

   findThis = 'AdjustSeq%/AdjCoilSensSeq';
   p = strfind(TextHeader.measYapsAscConv,findThis) + length(findThis);
   if Dim.Raw.nCha > 1 && ~isempty(p)
      Dim.Raw.nCha = Dim.Raw.nCha-1;
   end

   if Dim.Raw.nMaxRxChannels < Dim.Raw.nCha
      Dim.Raw.nCha = Dim.Raw.nMaxRxChannels;
   end

   if MrParcRaidFileCell{iMeas}.Control.File.turboFactor > 1
      Dim.Raw.nPhaseCor = 1;
   end
   if MrParcRaidFileCell{iMeas}.Control.File.epiFactor > 1
      Dim.Raw.nPhaseCor = 1;
   end
   %Control.User.readPhaseCor = Control.User.readPhaseCor && (Dim.Raw.nPhaseCor > 0);

   Dim.Recon.nNonOSColumns = round(Dim.Recon.nFourierColumns/Dim.Raw.readoutOversampleFactor);

   if Dim.Recon.nFourierPartitions == 1
      Dim.Recon.partitionFTLength = 1;
   end

   %% Raw Data Correction Factors (NOTE: CURRENTLY false!!!!)
   CoilSelectMap = containers.Map;
   findThis = '{\s*{\s*{\s*"[^"]+"[^\n]+';
   allCoilsInHeaderCell = regexp(TextHeader.complete, findThis, 'match');
   if ~isempty(allCoilsInHeaderCell)
      allCoilsInHeaderCell = allCoilsInHeaderCell{1};
      findThis = '{\s*{\s*"(?<name>[^"]+)"\s*}\s*{\s*(?<fft>[\d\.]+)\s*}\s*{\s*(?<re>[\d\.-]+)\s*}\s*{\s*(?<im>[\d\.-]+)\s*}\s*}';
      CoilStructArray = regexp(allCoilsInHeaderCell, findThis,'names');
      if length(CoilStructArray) == Dim.Raw.nCha
         for c=1:Dim.Raw.nCha
            tName = CoilStructArray(c).name;
            CoilSelect.fftScale = sscanf(CoilStructArray(c).fft,'%f');
            CoilSelect.rawDataCorrectionFactor = complex(sscanf(CoilStructArray(c).re,'%f'),sscanf(CoilStructArray(c).im,'%f'));
            CoilSelectMap(tName) = CoilSelect;
         end
      end
      if (length(keys(CoilSelectMap)) ~= Dim.Raw.nCha) && Control.User.applyRawDataCorrection
         if ~Control.User.silent
            display_message(Control.User.noGui,'ERROR','Non-unique channel names in CoilSelect','Error');
         end
         fclose(fid);
         Control = [];
         MrParcRaidFileCell = [];
         Dim = [];
         CoilSelectMap = [];
         CoilSelectMeasMap = [];
         RelativeSliceNumber = [];
         WiPMemBlock = [];
         return
      end
   end

   RelativeSliceNumber = [];
   if Control.User.readRelSliceNumber

      findThis = sprintf('ParamLong."relSliceNumber">[\\s\n]*{');
      pStart = regexp(TextHeader.complete,findThis);
      if ~isempty(pStart)
         pStart = pStart(1) + strfind(TextHeader.complete(pStart(1):pStart(1) + length(findThis)+20),'{');
         pStart = pStart(1);
         pEnd = strfind(TextHeader.complete(pStart(1):end),'}');
         if ~isempty(pEnd)
            pEnd = pStart + pEnd(1) - 2;
            findThis = '([^0-9\s]*[\-0-9]+)(\s+)';
            relSliceNumCell = regexp(TextHeader.complete(pStart:pEnd),findThis,'tokens');
            if ~isempty(relSliceNumCell)
               relativeSliceNumber = zeros(1,length(relSliceNumCell));
               for i=1:length(relSliceNumCell)
                  relativeSliceNumber(i) = str2double(relSliceNumCell{i}{1});
               end
               RelativeSliceNumber{iMeas} = relativeSliceNumber;
               Control.User.foundRelSliceNumber = true;
            end
         end
      end
   end

   % WiPMemBlock
   findThis = sprintf('ParamMap."sWiPMemBlock">[\\s\n]*{');
   pStartWiPMemBlock = regexp(TextHeader.complete, findThis, 'ignorecase');
   WiPMemBlock.adRes = [];
   WiPMemBlock.tFree = [];
   WiPMemBlock.adFree = [];
   WiPMemBlock.lCSatBW = [];
   WiPMemBlock.alFree = [];
   if ~isempty(pStartWiPMemBlock)
      pStartWiPMemBlock = pStartWiPMemBlock(1) + strfind(TextHeader.complete(pStartWiPMemBlock(1):pStartWiPMemBlock(1) + length(findThis)+20),'{');
      pStartWiPMemBlock = pStartWiPMemBlock(1);
      findThis = sprintf('<ParamLong."alFree">[\\s\n]*{');
      pStartAlFree = regexp(TextHeader.complete(pStartWiPMemBlock:end), findThis, 'ignorecase');
      if ~isempty(pStartAlFree)
         pStartAlFree = pStartWiPMemBlock + pStartAlFree(1); 
         pStartAlFree = pStartAlFree + strfind(TextHeader.complete(pStartAlFree(1):pStartAlFree(1) + length(findThis)+20),'{');
         pStartAlFree = pStartAlFree(1);
         pEnd_alFree = strfind(TextHeader.complete(pStartAlFree(1):end),'}');
         if ~isempty(pEnd_alFree)
            pEnd_alFree = pStartAlFree + pEnd_alFree(1) - 2;
            % Skip any precision/default stuff
            pLessThan = strfind(TextHeader.complete(pStartAlFree:pEnd_alFree),'<');
            if ~isempty(pLessThan)
               pStartAlFree = pLessThan(end) + pStartAlFree;
               pStartAlFree = pStartAlFree + strfind(TextHeader.complete(pStartAlFree:pEnd_alFree),sprintf('\n'));
               pStartAlFree = pStartAlFree(1);
            end
            alFreeCell = regexp(TextHeader.complete(pStartAlFree:pEnd_alFree),'([\-0-9]+)\s+','tokens');
            if ~isempty(alFreeCell)
               Control.User.foundWiPMemBlock = true;
               WiPMemBlock.alFree = zeros(1,length(alFreeCell),'int32');
               for i_alFree=1:length(alFreeCell)
                  WiPMemBlock.alFree(i_alFree) = int32(str2num(alFreeCell{i_alFree}{1}));
               end
            end
         end
      end
      findThis = sprintf('<ParamLong."lCSatBW">[\\s\n]*{');
      pStartLCSatBW = regexp(TextHeader.complete(pStartWiPMemBlock:end), findThis, 'ignorecase');
      if ~isempty(pStartLCSatBW)
         pStartLCSatBW = pStartWiPMemBlock + pStartLCSatBW(1); 
         pStartLCSatBW = pStartLCSatBW + strfind(TextHeader.complete(pStartLCSatBW(1):pStartLCSatBW(1) + length(findThis)+20),'{');
         pStartLCSatBW = pStartLCSatBW(1);
         pEnd_lCSatBW = strfind(TextHeader.complete(pStartLCSatBW(1):end),'}');
         if ~isempty(pEnd_lCSatBW)
            pEnd_lCSatBW = pStartLCSatBW + pEnd_lCSatBW(1) - 2;
            % Skip any precision/default stuff
            pLessThan = strfind(TextHeader.complete(pStartLCSatBW:pEnd_lCSatBW),'<');
            if ~isempty(pLessThan)
               pStartLCSatBW = pLessThan(end) + pStartLCSatBW;
               pStartLCSatBW = pStartLCSatBW + strfind(TextHeader.complete(pStartLCSatBW:pEnd_lCSatBW),sprintf('\n'));
               pStartLCSatBW = pStartLCSatBW(1);
            end
            lCSatBWCell = regexp(TextHeader.complete(pStartLCSatBW:pEnd_lCSatBW),'([\-0-9]+)\s+','tokens');
            if ~isempty(lCSatBWCell)
               WiPMemBlock.lCSatBW = zeros(1,length(lCSatBWCell),'int32');
               Control.User.foundWiPMemBlock = true;
               for i_lCSatBW=1:length(lCSatBWCell)
                  WiPMemBlock.lCSatBW(i_lCSatBW) = int32(str2num(lCSatBWCell{i_lCSatBW}{1}));
               end
            end
         end
      end
      findThis = sprintf('<ParamDouble."adFree">[\\s\n]*{');
      pStartAdFree = regexp(TextHeader.complete(pStartWiPMemBlock:end), findThis, 'ignorecase');
      if ~isempty(pStartAdFree)
         pStartAdFree = pStartWiPMemBlock + pStartAdFree(1); 
         pStartAdFree = pStartAdFree + strfind(TextHeader.complete(pStartAdFree(1):pStartAdFree(1) + length(findThis)+20),'{');
         pStartAdFree = pStartAdFree(1);
         pEnd_adFree = strfind(TextHeader.complete(pStartAdFree(1):end),'}');
         if ~isempty(pEnd_adFree)
            pEnd_adFree = pStartAdFree + pEnd_adFree(1) - 2;
            % Skip any precision/default stuff
            pLessThan = strfind(TextHeader.complete(pStartAdFree:pEnd_adFree),'<');
            if ~isempty(pLessThan)
               pStartAdFree = pLessThan(end) + pStartAdFree;
               pStartAdFree = pStartAdFree + strfind(TextHeader.complete(pStartAdFree:pEnd_adFree),sprintf('\n'));
               pStartAdFree = pStartAdFree(1);
            end
            adFreeCell = regexp(TextHeader.complete(pStartAdFree:pEnd_adFree),'([\-0-9\.]+)\s+','tokens');
            if ~isempty(adFreeCell)
               Control.User.foundWiPMemBlock = true;
               WiPMemBlock.adFree = zeros(1,length(adFreeCell));
               for i_adFree=1:length(adFreeCell)
                  WiPMemBlock.adFree(i_adFree) = str2num(adFreeCell{i_adFree}{1});
               end
            end
         end
      end
      findThis = sprintf('<ParamDouble."adRes">[\\s\n]*{');
      pStartAdRes = regexp(TextHeader.complete(pStartWiPMemBlock:end), findThis, 'ignorecase');
      if ~isempty(pStartAdRes)
         pStartAdRes = pStartWiPMemBlock + pStartAdRes(1); 
         pStartAdRes = pStartAdRes + strfind(TextHeader.complete(pStartAdRes(1):pStartAdRes(1) + length(findThis)+20),'{');
         pStartAdRes = pStartAdRes(1);
         pEnd_adRes = strfind(TextHeader.complete(pStartAdRes(1):end),'}');
         if ~isempty(pEnd_adRes)
            pEnd_adRes = pStartAdRes + pEnd_adRes(1) - 2;
            % Skip any precision/default stuff
            pLessThan = strfind(TextHeader.complete(pStartAdRes:pEnd_adRes),'<');
            if ~isempty(pLessThan)
               pStartAdRes = pLessThan(end) + pStartAdRes;
               pStartAdRes = pStartAdRes + strfind(TextHeader.complete(pStartAdRes:pEnd_adRes),sprintf('\n'));
               pStartAdRes = pStartAdRes(1);
            end
            adResCell = regexp(TextHeader.complete(pStartAdRes:pEnd_adRes),'([\-0-9\.]+)\s+','tokens');
            if ~isempty(adResCell)
               Control.User.foundWiPMemBlock = true;
               WiPMemBlock.adRes = zeros(1,length(adResCell));
               for i_adRes=1:length(adResCell)
                  WiPMemBlock.adRes(i_adRes) = str2num(adResCell{i_adRes}{1});
               end
            end
         end
      end
      findThis = sprintf('<ParamString."tFree">[\\s\n]*{');
      pStartTFree = regexp(TextHeader.complete(pStartWiPMemBlock:end), findThis, 'ignorecase');
      if ~isempty(pStartTFree)
         pStartTFree = pStartWiPMemBlock + pStartTFree(1); 
         pStartTFree = pStartTFree + strfind(TextHeader.complete(pStartTFree(1):pStartTFree(1) + length(findThis)+20),'{');
         pStartTFree = pStartTFree(1);
         pEnd_tFree = strfind(TextHeader.complete(pStartTFree(1):end),'}');
         if ~isempty(pEnd_tFree)
            pEnd_tFree = pStartTFree + pEnd_tFree(1) - 2;
            % Skip any precision/default stuff
            pLessThan = strfind(TextHeader.complete(pStartTFree:pEnd_tFree),'<');
            if ~isempty(pLessThan)
               pStartTFree = pLessThan(end) + pStartTFree;
               pStartTFree = pStartTFree + strfind(TextHeader.complete(pStartTFree:pEnd_tFree),sprintf('\n'));
               pStartTFree = pStartTFree(1);
            end
            tFreeCell = regexp(TextHeader.complete(pStartTFree:pEnd_tFree),'"([^"])"\s+','tokens');
            if ~isempty(tFreeCell)
               Control.User.foundWiPMemBlock = true;
               WiPMemBlock.tFree = zeros(1,length(tFreeCell));
               for i_tFree=1:length(tFreeCell)
                  WiPMemBlock.tFree(i_tFree) = str2num(tFreeCell{i_tFree}{1});
               end
            end
         end
      end
   end % if wipmemblock found

   %% FFT Correction Factors
   for c=1:Dim.Raw.nCha
      findThis = 'asList\[';
      findThis = sprintf('%s%d',findThis,c-1);
      findThis = [findThis,'\]\.sCoilElementID.tElement\s*=\s*"([^"]+)"\s*'];
      p = regexp(TextHeader.measYapsAscConv,findThis,'tokens');
      CoilSelectMeas.tElement = p{1}{1};
      findThis = 'asList\[';
      findThis = sprintf('%s%d',findThis,c-1);
      findThis = [findThis,'\]\.lRxChannelConnected\s*=\s*([0-9.]+)\s*'];
      p = regexp(TextHeader.measYapsAscConv,findThis,'tokens');
      CoilSelectMeas.lRxChannel = str2num(p{1}{1});
      findThis = 'asList\[';
      findThis = sprintf('%s%d',findThis,c-1);
      findThis = [findThis,'\]\.lADCChannelConnected\s*=\s*([0-9.]+)\s*'];
      p = regexp(TextHeader.measYapsAscConv,findThis,'tokens');
      lADCChannel = str2num(p{1}{1});
      if MrParcRaidFileCell{iMeas}.Control.User.resetFFTScale
         CoilSelectMeas.flFFTCorrectionFactor = single(1.0); 
      else
         findThis = 'aFFT_SCALE\[';
         findThis = sprintf('%s%d',findThis,c-1);
         findThis = [findThis,'\]\.flFactor\s*=\s*([0-9.]+)\s*'];
         p = regexp(TextHeader.measYapsAscConv,findThis,'tokens');
         CoilSelectMeas.flFFTCorrectionFactor = single(str2double(p{1}{1}));
      end
      keySet(c) = lADCChannel;
      valueSet(c) = CoilSelectMeas;
   end
   if Dim.Raw.nCha > 1
      CoilSelectMeasMap = containers.Map(keySet,arrayfun(@(x) ({x}), valueSet));
   else
      CoilSelectMeasMap = containers.Map(keySet,valueSet);
   end

   clear keySet valueSet TextHeader.measYapsAscConv;

   if Dim.Recon.nFourierLines == Dim.Recon.nLin && Dim.Recon.nFourierPartitions == Dim.Recon.nPar
      MrParcRaidFileCell{iMeas}.Control.User.noFillToFullFourier = true;
   end

   if ~MrParcRaidFileCell{iMeas}.Control.File.is3D
      Dim.Recon.nPar = Dim.Recon.nFourierPartitions;
   end

   clear TextHeader;

   if Dim.Recon.nFourierLines > Dim.Recon.nLin
       Dim.Recon.nLin = Dim.Recon.nFourierLines;
   end

   if Dim.Recon.nFourierPartitions > Dim.Recon.nPar
       Dim.Recon.nPar = Dim.Recon.nFourierPartitions;
   end

end

% =============================================================================
function textMessage = print_setup_status(filenameIn, UserInputDefault, Control)
% -----------------------------------------------------------------------------
% Send text to the user based on GUI or print statement.
% -----------------------------------------------------------------------------
   textMessage = sprintf('Filename=%s', filenameIn);
   for iArg = 1:length(UserInputDefault)
      textMessage = sprintf('%s\n%s=%d', ...
                            textMessage, ...
                            UserInputDefault(iArg).varArgIn, ...
                            Control.User.(UserInputDefault(iArg).fieldName));
   end
   display_message(Control.User.noGui,'WAIT',textMessage,'File Will Be Read With These Parameters');
end

% =============================================================================
function textMessage = print_available_options(UserInputDefault)
% -----------------------------------------------------------------------------
% Based on default inputs, indicate to user what is possible.
% -----------------------------------------------------------------------------
   textMessage = sprintf('[kSpace, Other] = fast_read_vd13(fileName, ''option1'', ''option2'', ...');
   textMessage = sprintf('%s\n\nfileName is the name of a meas.dat to read and',textMessage);
   textMessage = sprintf('%s\n\nOptions are as follows:', textMessage);
   for iArg = 1:length(UserInputDefault)
      textMessage = sprintf('%s\n\n''%s'': (default', ...
                            textMessage, ...
                            UserInputDefault(iArg).varArgIn);
      switch UserInputDefault(iArg).type
         case 'logical'
            if UserInputDefault(iArg).defaultValue
               textMessage = sprintf('%s true)', textMessage);
            else
               textMessage = sprintf('%s false)', textMessage);
            end
         case 'uint64'
            textMessage = sprintf('%s %d)', textMessage, UserInputDefault(iArg).defaultValue);
         otherwise
            textMessage = sprintf('%s unknown)');
      end
      textMessage = sprintf('%s %s', textMessage, UserInputDefault(iArg).description);
   end
   if usejava('jvm') && ~feature('ShowFigureWindows')
      display(textMessage);
   else
      msgbox(textMessage, 'Usage');
   end
end

% =============================================================================
function MdhBlockMeta = append_to_mdhblockmeta(iBlock, iBlockMeta, MdhBlock, varargin)
% -----------------------------------------------------------------------------
% Append to existing MdhBlockMeta list.
% -----------------------------------------------------------------------------
   copyField = {'isSyncData', 'isAcqEnd', 'length', 'nScan', 'nSamplesInScan', 'nChannelUsed'};
   if length(varargin) > 0
      MdhBlockMeta = varargin{1};
   end
   MdhBlockMeta(iBlockMeta).nBlock = 1;
   MdhBlockMeta(iBlockMeta).iBlock(MdhBlockMeta(iBlockMeta).nBlock) = iBlock;
   for iField=1:length(copyField)
      MdhBlockMeta(iBlockMeta).(copyField{iField}) = MdhBlock(iBlock).(copyField{iField});
   end
end

% =============================================================================
function set_global_mdh_parameters(filenameIn)
% -----------------------------------------------------------------------------
% Handle some global parameters.
% -----------------------------------------------------------------------------

   global MDH_SCANSIZE MDH_CHANSIZE MDH_PMUSIZE MAX_BYTES_PER_FREAD
   global mdhColumns mdhBitFlag

   % Constants
   MDH_ACQEND = 2^uint32(0);
   MDH_RTFEEDBACK = 2^uint32(1);
   MDH_SYNCDATA = 2^uint32(5);
   MDH_RAWDATACORRECTION = 2^uint32(10);
   MDH_PHASCOR= 2^uint32(21);
   MDH_NOISEADJSCAN = 2^uint32(25);
   MDH_REFLECT = 2^uint32(24);

   % Some default sizes
   MDH_SCANSIZE = uint64(192);
   MDH_CHANSIZE = uint64(32);
   MDH_PMUSIZE = uint64(60);
   MAX_BYTES_PER_FREAD = 2^uint64(27);
   fileSpecs = dir(filenameIn);
   MAX_BYTES_PER_FREAD = min(MAX_BYTES_PER_FREAD,fileSpecs.bytes); % Constrain to file size

   % Columns within the binary scan lines where
   %    particular data is found along with how
   %    it should be converted.
   mdhColumns.scan.dmaLength.c = [1 4];
   mdhColumns.scan.dmaLength.t = 'uint32';
   mdhColumns.scan.timeStamp.c = [13 16];
   mdhColumns.scan.timeStamp.t = 'uint32';
   mdhColumns.scan.evalMask1.c = [41 44];
   mdhColumns.scan.evalMask1.t = 'uint32';
   mdhColumns.scan.evalMask2.c = [45 48];
   mdhColumns.scan.evalMask2.t = 'uint32';
   mdhColumns.scan.nSamplesInScan.c = [49 50];
   mdhColumns.scan.nSamplesInScan.t = 'uint16';
   mdhColumns.scan.nChannelUsed.c = [51 52];
   mdhColumns.scan.nChannelUsed.t = 'uint16';
   mdhColumns.scan.indexLIN.c = [53 54];
   mdhColumns.scan.indexLIN.t = 'uint16';
   mdhColumns.scan.indexACQ.c = [55 56];
   mdhColumns.scan.indexACQ.t = 'uint16';
   mdhColumns.scan.indexSLC.c = [57 58];
   mdhColumns.scan.indexSLC.t = 'uint16';
   mdhColumns.scan.indexPAR.c = [59 60];
   mdhColumns.scan.indexPAR.t = 'uint16';
   mdhColumns.scan.indexECO.c = [61 62];
   mdhColumns.scan.indexECO.t = 'uint16';
   mdhColumns.scan.indexPHS.c = [63 64];
   mdhColumns.scan.indexPHS.t = 'uint16';
   mdhColumns.scan.indexREP.c = [65 66];
   mdhColumns.scan.indexREP.t = 'uint16';
   mdhColumns.scan.indexSET.c = [67 68];
   mdhColumns.scan.indexSET.t = 'uint16';
   mdhColumns.scan.indexSEG.c = [69 70];
   mdhColumns.scan.indexSEG.t = 'uint16';
   mdhColumns.scan.cutOffDataPre.c = [81 82];
   mdhColumns.scan.cutOffDataPre.t = 'uint16';
   mdhColumns.scan.cutOffDataPost.c = [83 84];
   mdhColumns.scan.cutOffDataPost.t = 'uint16';
   mdhColumns.scan.kSpaceCentreLineNo.c = [97 98];
   mdhColumns.scan.kSpaceCentreLineNo.t = 'uint16';
   mdhColumns.scan.kSpaceCentrePartitionNo.c = [99 100];
   mdhColumns.scan.kSpaceCentrePartitionNo.t = 'uint16';

   mdhColumns.chan.sizeChanInBytes.c = [1 4];
   mdhColumns.chan.sizeChanInBytes.t = 'uint32';
   mdhColumns.chan.channelID.c = [25 26];
   mdhColumns.chan.channelID.t = 'uint16';

   mdhColumns.pmu.timeStamp.c = [1 8];
   mdhColumns.pmu.timeStamp.t = 'uint32';
   mdhColumns.pmu.counter.c = [9 12];
   mdhColumns.pmu.counter.t = 'uint32';
   mdhColumns.pmu.duration.c = [13 16];
   mdhColumns.pmu.duration.t = 'uint32';
   mdhColumns.pmu.data.c = [17 0];
   mdhColumns.pmu.data.t = 'uint32';

   mdhBitFlag.isNoiseAdj.s = 'evalMask1';
   mdhBitFlag.isNoiseAdj.f = MDH_NOISEADJSCAN;
   mdhBitFlag.isPhaseCor.s = 'evalMask1';
   mdhBitFlag.isPhaseCor.f = MDH_PHASCOR;
   mdhBitFlag.isRTFeedback.s = 'evalMask1';
   mdhBitFlag.isRTFeedback.f = MDH_RTFEEDBACK;
   mdhBitFlag.isAcqEnd.s = 'evalMask1';
   mdhBitFlag.isAcqEnd.f = MDH_ACQEND;
   mdhBitFlag.isReflect.s = 'evalMask1';
   mdhBitFlag.isReflect.f = MDH_REFLECT;
   mdhBitFlag.isRawDataCorrection.s = 'evalMask1';
   mdhBitFlag.isRawDataCorrection.f = MDH_RAWDATACORRECTION;
   mdhBitFlag.isSyncData.s = 'evalMask1';
   mdhBitFlag.isSyncData.f = MDH_SYNCDATA;

end

% =============================================================================
function MdhParam = extract_single_mdh_parameter(mdhBuffer, paramType, nameToExtract)
% -----------------------------------------------------------------------------
% Extract a parameter from the uint8 array based on input description.
% -----------------------------------------------------------------------------

   global mdhColumns mdhBitFlag

   columns = mdhColumns.(paramType).(nameToExtract).c;
   temp = mdhBuffer(columns(1):columns(2),:);
   MdhParam = typecast(temp(:),mdhColumns.(paramType).(nameToExtract).t);

   return

end

% =============================================================================
function MdhParam = extract_mdh_parameters(mdhBuffer, paramType, namesToExtract)
% -----------------------------------------------------------------------------
% Extract an array of parameters, similar to extract_single_mdh_parameter.
% -----------------------------------------------------------------------------

   global mdhColumns mdhBitFlag

   switch nargin
      case 3
         fieldNames = namesToExtract;
      case 2
         fieldNames = fieldnames(mdhColumns.(paramType));
      otherwise
         paramType = 'scan';
         fieldNames = fieldnames(mdhColumns.(paramType));
   end

   for iField = 1:length(fieldNames)
      columns = mdhColumns.(paramType).(fieldNames{iField}).c;
      if columns(2) == 0 % PMU data goes to the end of the buffer
         columns(2) = length(mdhBuffer);
      end
      temp = mdhBuffer(columns(1):columns(2),:);
      MdhParam.(fieldNames{iField}) = typecast(temp(:),mdhColumns.(paramType).(fieldNames{iField}).t);
   end

   return

end

% =============================================================================
function MdhBitFlags = extract_mdh_bitflags(mdhBuffer, flagsToExtract)
% -----------------------------------------------------------------------------
% Pull bitflags from the MDH evalMasks
% -----------------------------------------------------------------------------

   global mdhColumns mdhBitFlag
   
   if nargin == 1
      fieldNames = fieldnames(mdhBitFlag);
   else
      fieldNames = flagsToExtract;
   end

   masks = {'evalMask1', 'evalMask2'};
   for iMask = 1:2
      MdhBitFlags.(masks{iMask}) = extract_single_mdh_parameter(mdhBuffer, 'scan', masks{iMask});
   end
   for iField = 1:length(fieldNames)
      evalMask = mdhBitFlag.(fieldNames{iField}).s;
      bitFlag = mdhBitFlag.(fieldNames{iField}).f;
      MdhBitFlags.(fieldNames{iField}) = bitand(MdhBitFlags.(evalMask), bitFlag) == bitFlag;
   end

   return;

end

% =============================================================================
function MdhBlock = init_mdh_block(offsetInFile, dmaLength, MdhBitFlags, MdhParam)
% -----------------------------------------------------------------------------
% Set up MdhBlock structure
% -----------------------------------------------------------------------------

   MdhBlock.offsetInFile = offsetInFile;
   MdhBlock.length = dmaLength;
   MdhBlock.nScan = uint64(0); % Initialize with zero value
   MdhBlock.isSyncData = MdhBitFlags.isSyncData;
   MdhBlock.isAcqEnd = MdhBitFlags.isAcqEnd;
   MdhBlock.nSamplesInScan = uint64(MdhParam.nSamplesInScan);
   MdhBlock.nChannelUsed = MdhParam.nChannelUsed;

end

% =============================================================================
function dimOut = max_from_mdh_param(dimIn, nChannelUsed, dimNames, MdhParam, indices)
% -----------------------------------------------------------------------------
% Compute output dimensions
% -----------------------------------------------------------------------------
   dimOut = dimIn;
   dimOut(1) = max(dimIn(1), nChannelUsed-1);
   for iDim=2:numel(dimNames)
      indexName = strcat('index',dimNames{iDim});
      dimOut(iDim) = max(dimIn(iDim), max(MdhParam.(indexName)(indices(:))));
   end
   return
end

% =============================================================================
function rankOut = rank_from_mdh_param(rankIn, dimNames, MdhParam, nChannels, indices)
% -----------------------------------------------------------------------------
% Determine which dimensions vary most rapidly based on mdh values
% -----------------------------------------------------------------------------
   rankOut = rankIn;
   rankOut(1) = rankIn(1) + numel(indices(2:end)) * nChannels;
   for iDim=2:numel(dimNames)
      indexName = strcat('index',dimNames{iDim});
      diffMask = MdhParam.(indexName)(indices(2:end)) ~= MdhParam.(indexName)(indices(1:end-1));
      rankOut(iDim) = rankIn(iDim) + sum(diffMask,1);
   end
   return
end

% =============================================================================
function rankOut = rank_from_index_prev(rankIn, dimNames, MdhParam, nChannels, index, indexPrev)
% -----------------------------------------------------------------------------
% Determine column order based on indices in the mdh
% -----------------------------------------------------------------------------
   rankOut = rankIn;
   rankOut(1) = rankOut(1) + nChannels;
   for iDim=2:numel(dimNames)
      indexName = strcat('index',dimNames{iDim});
      diffMask = indexPrev.(indexName) ~= MdhParam.(indexName)(index);
      rankOut(iDim) = rankIn(iDim) + sum(diffMask,1);
   end
   return
end

% =============================================================================
function indexPrev = set_index_prev(dimNames, MdhParam, index)
% -----------------------------------------------------------------------------
% Bookeeping for tracking what dimensions vary most rapidly
% -----------------------------------------------------------------------------
   for iDim=2:numel(dimNames)
      indexName = strcat('index',dimNames{iDim});
      indexPrev.(indexName) = MdhParam.(indexName)(index);
   end
   return
end
 
% =============================================================================
function result = is_field_in_struct(structIn, fieldName)
% -----------------------------------------------------------------------------
% Determine if the input field name is a leaf/field in a struct
% -----------------------------------------------------------------------------
   if ~isstruct(structIn)
      result = 0;
      return
   end
   result = sum(strcmp(fieldnames(structIn), fieldName)) > 0;
   return
end


% =============================================================================
function output = compute_cumulative_product(input)
% -----------------------------------------------------------------------------
% Compute the cumulative product for an input array.  Apparently before matlab
%    2013, cumprod only supports floating cumprod and I need one for integers!
%
% output = compute_cumulative_product(input) returns array of same type and
%    structure as input with cumulative product computed from first to last
%    element
%
% HISTORY:   05 Dec 2014. J. Roberts. roberts@ucair.med.utah.edu
%            UCAIR, Dept Radiology, School of Med, U of U, SLC, UT
% -----------------------------------------------------------------------------

   output = ones(size(input),class(input(1)));
   output(1) = input(1);
   for p=2:numel(input)
      output(p) = output(p-1) * input(p);
   end

end

% =============================================================================
function output = compute_product(input)
% -----------------------------------------------------------------------------
% Compute the product for an input array.  Apparently before matlab 2013, prod
%    only supports floating prod and I need one for integers!
%
% output = compute_product(input) returns array of same type and structure
%                                     as input with product 
%                                     computed from first to last element
%
% HISTORY:   05 Dec 2014. J. Roberts. roberts@ucair.med.utah.edu
%            UCAIR, Dept Radiology, School of Med, U of U, SLC, UT
% -----------------------------------------------------------------------------

   output = input(1);
   for p=2:numel(input)
      output = output * input(p);
   end

end

% =============================================================================
function display_message(noGui,typeMessage,textMessage,textExtra)
% -----------------------------------------------------------------------------

   if noGui
      disp(textMessage);
   else
      if strcmp(typeMessage,'WARNING')
         warndlg(textMessage);
      elseif strcmp(typeMessage,'ERROR')
         errordlg(textMessage,textExtra);
      else
         uiwait(msgbox(textMessage, textExtra, 'modal'));
      end
   end

end

% =============================================================================
function [sizeAllScansInBytes, sizeChanInBytes, MdhBlock, nBlockSYNCDATA, ...
          nBlockACQEND, nBlockScan, blockNScan, nScanInMeas, dmaLengthScan, ...
          nSamplesInScan] = chunk_into_equal_line_length_blocks(fid,MrParcRaidFileCell, iMeas, Dim, Control)
% -----------------------------------------------------------------------------
% -----------------------------------------------------------------------------
   global MDH_SCANSIZE MDH_CHANSIZE MDH_PMUSIZE MAX_BYTES_PER_FREAD

   sizeAllScansInBytes = MrParcRaidFileCell{iMeas}.totalSize - Dim.ullMeasHeaderSize;
   fseek(fid,MrParcRaidFileCell{iMeas}.offsetInFile + Dim.ullMeasHeaderSize,'bof');
   % Read the first line to get the DMA length and line type
   mdhScan = fread(fid,MDH_SCANSIZE,'uchar=>uchar');
   MdhParam = extract_mdh_parameters(mdhScan, 'scan', {'dmaLength','nSamplesInScan','nChannelUsed'});
   MdhBitFlags = extract_mdh_bitflags(mdhScan, {'isAcqEnd', 'isSyncData'});
   bHasACQEND = MdhBitFlags.isAcqEnd;
   dmaLength = uint64(bitand(MdhParam.dmaLength,hex2dec('FFFFFF'))); % Only first 3 bytes used for dmaLength
   iRelativePosition = uint64(1);
   iBlock = 1;
   offsetInFile = MrParcRaidFileCell{iMeas}.offsetInFile + Dim.ullMeasHeaderSize;
   % NOTE: We make a very simple assumption that the same number of
   %    samples per scan will hold for scans of the same dmaLength
   % Similarly for the number of channels (which mean nothing in SYNCDATA scans)
   MdhBlock(iBlock) = init_mdh_block(offsetInFile, dmaLength, MdhBitFlags, MdhParam);
   while ((iRelativePosition < sizeAllScansInBytes) && ~bHasACQEND)
      if ~(MrParcRaidFileCell{iMeas}.Control.User.noGui || MrParcRaidFileCell{iMeas}.Control.User.silent)
         textMessage = 'Pass 1 of 3: Parsing blocks';
         waitbar(double(iRelativePosition)/double(sizeAllScansInBytes), Control.GUI.hWaitBar, textMessage);
      end
      fseek(fid,MdhBlock(iBlock).offsetInFile + MdhBlock(iBlock).nScan * dmaLength,'bof');
      bytesRemaining = sizeAllScansInBytes-iRelativePosition+uint64(1);
      if dmaLength < 2048
         nScanToRead=uint64(10);  % We don'p usually get a lot short lines
      else
         nScanToRead = max(idivide(MAX_BYTES_PER_FREAD, dmaLength,'fix'),uint64(1));
      end   
      if (nScanToRead*dmaLength > bytesRemaining)
         nScanToRead = idivide(bytesRemaining,dmaLength,'fix');
         if nScanToRead < 1
            mdhScan = fread(fid,MDH_SCANSIZE,'uchar=>uchar');
            MdhParam = extract_mdh_parameters(mdhScan, 'scan', {'dmaLength','nSamplesInScan','nChannelUsed'});
            MdhBitFlags = extract_mdh_bitflags(mdhScan, {'isAcqEnd', 'isSyncData'});
            dmaLengthNext = uint64(bitand(MdhParam.dmaLength,hex2dec('FFFFFF'))); % Only first 3 bytes used by dmaLength
            nScanToRead = uint64(1);
            if dmaLengthNext == dmaLength
               error_message
            else
               iBlock = iBlock + 1;
               offsetInFile = MdhBlock(iBlock-1).offsetInFile + MdhBlock(iBlock-1).nScan * dmaLength;
               dmaLength = dmaLengthNext;
               MdhBlock(iBlock) = init_mdh_block(offsetInFile, dmaLength, MdhBitFlags, MdhParam);
               fseek(fid,MdhBlock(iBlock).offsetInFile + MdhBlock(iBlock).nScan * dmaLength,'bof');
            end
         end
      end
      [blockOfScans readCount] = fread(fid,dmaLength*nScanToRead,'uchar=>uchar');
      if MdhBlock(iBlock).isAcqEnd && (readCount < (dmaLength*nScanToRead))
         dmaLength = readCount / nScanToRead;
      end
      for iStart=1:dmaLength:dmaLength*nScanToRead
         mdhScan = blockOfScans(iStart:iStart+dmaLength-1);
         MdhParam = extract_mdh_parameters(mdhScan, 'scan', {'dmaLength','nSamplesInScan','nChannelUsed'});
         MdhBitFlags = extract_mdh_bitflags(mdhScan, {'isAcqEnd', 'isSyncData'});
         dmaLengthNext = uint64(bitand(MdhParam.dmaLength,hex2dec('FFFFFF'))); % Only first 3 bytes for dmaLength
         bHasACQEND = bHasACQEND || MdhBlock(iBlock).isAcqEnd;
         if dmaLengthNext == dmaLength
            MdhBlock(iBlock).nScan = MdhBlock(iBlock).nScan + 1;
            iRelativePosition = iRelativePosition + dmaLength;
         else
            iBlock = iBlock + 1;
            offsetInFile =MdhBlock(iBlock-1).offsetInFile + MdhBlock(iBlock-1).nScan * dmaLength;
            dmaLength = dmaLengthNext;
            MdhBlock(iBlock) = init_mdh_block(offsetInFile, dmaLength, MdhBitFlags, MdhParam);
            break;
         end
         if bHasACQEND
            break;
         end
      end % for over block of scans with same dmalength
   end % while over full file
   if ~(MrParcRaidFileCell{iMeas}.Control.User.noGui || MrParcRaidFileCell{iMeas}.Control.User.silent)
      textMessage = 'Pass 1 of 3: Parsing blocks';
      waitbar(double(iRelativePosition)/double(sizeAllScansInBytes), Control.GUI.hWaitBar, textMessage);
   end
   if ((iRelativePosition < sizeAllScansInBytes) && bHasACQEND)         
      if ~MrParcRaidFileCell{iMeas}.Control.User.silent
         display_message(MrParcRaidFileCell{iMeas}.Control.User.noGui,'WARNING','ACQ_END before end of file.');
      end
   end
   if ~bHasACQEND
      if ~MrParcRaidFileCell{iMeas}.Control.User.silent
         display_message(MrParcRaidFileCell{iMeas}.Control.User.noGui,'WARNING','No ACQ_END encountered.  Data possibly incomplete or corrupted.');
      end
   end
   nBlockSYNCDATA = 0;
   nBlockACQEND = 0;
   nBlockScan=0;
   for iBlock=1:length(MdhBlock)
      if MdhBlock(iBlock).isAcqEnd
         nBlockACQEND = nBlockACQEND + 1;
      elseif MdhBlock(iBlock).isSyncData
         nBlockSYNCDATA = nBlockSYNCDATA + 1;
      else
         nBlockScan = nBlockScan + 1;
         blockNScan(nBlockScan) = MdhBlock(iBlock).nScan;
         dmaLengthScan(nBlockScan) = MdhBlock(iBlock).length;
         nSamplesInScan(nBlockScan) = MdhBlock(iBlock).nSamplesInScan;
      end
   end
   if nBlockScan < 1
      if ~MrParcRaidFileCell{iMeas}.Control.User.silent
         display_message(MrParcRaidFileCell{iMeas}.Control.User.noGui,'ERROR','No readouts detected.','Error');
      end
      fclose(fid);
      sizeAllScansInBytes = [];
      sizeChanInBytes = [];
      MdhBlock = [];
      nBlockSYNCDATA = [];
      nBlockACQEND = [];
      nBlockScan = [];
      blockNScan = [];
      nScanInMeas = [];
      dmaLengthScan = [];
      nSamplesInScan = [];
      return
   end
   nScanInMeas = uint64(sum(blockNScan(:)));
   if (length(unique(nSamplesInScan)) ~= 1) || (length(unique(dmaLengthScan)) ~= 1)
      if ~MrParcRaidFileCell{iMeas}.Control.User.silent
         display_message(MrParcRaidFileCell{iMeas}.Control.User.noGui,'ERROR','Multiple different samples per scan detected.','Error');
      end
      fclose(fid);
      sizeAllScansInBytes = [];
      sizeChanInBytes = [];
      MdhBlock = [];
      nBlockSYNCDATA = [];
      nBlockACQEND = [];
      nBlockScan = [];
      blockNScan = [];
      nScanInMeas = [];
      dmaLengthScan = [];
      nSamplesInScan = [];
      return
   end
   nSamplesInScan = unique(nSamplesInScan);
   nSamplesInScan = nSamplesInScan(1);
   % Read the first real scan channel header to get scan length in bytes
   for iBlock=1:length(MdhBlock)
      if MdhBlock(iBlock).isAcqEnd || MdhBlock(iBlock).isSyncData
         continue
      else
         fseek(fid,MdhBlock(iBlock).offsetInFile+MDH_SCANSIZE,'bof');
         mdhChan = fread(fid,MDH_CHANSIZE,'uchar=>uchar');
         sizeChanInBytes = uint64(bitshift(extract_single_mdh_parameter(mdhChan,'chan','sizeChanInBytes'),-8));
         break
      end
   end

end

% =============================================================================
function PMUOut = extract_pmu(fid,MdhBlock)
% -----------------------------------------------------------------------------
% Read the PMU data from SYNCDATA lines
% -----------------------------------------------------------------------------

   global MDH_SCANSIZE MDH_CHANSIZE MDH_PMUSIZE MAX_BYTES_PER_FREAD

   PMU.TimeStamp.nColumn = 0;
   PMU.TimeStamp.raw = reshape([],2,0);
   PMU.TimeStamp.data = reshape([],2,0);
   Waveform.typeName = '';
   Waveform.nColumn = 0;
   Waveform.nRow = 0;
   Waveform.period = 0;
   Waveform.data = [];
   Waveform.trigger = [];
   pmuTypes = {'ECG1', 'ECG2', 'ECG3', 'ECG4', 'PULS', 'RESP', 'EXT1', 'EXT2'};
   PMU.Waveforms = repmat(Waveform,1,length(pmuTypes));
   for iType=1:length(pmuTypes)
      PMU.Waveforms(iType).typeName = pmuTypes(iType);
   end
   for iBlock = 1:length(MdhBlock)

      % Skip non-SYNCDATA
      if ~MdhBlock(iBlock).isSyncData 
         continue
      end

      fseek(fid,MdhBlock(iBlock).offsetInFile,'bof');
      nScanToRead = MdhBlock(iBlock).nScan;
      blockOfScans = fread(fid,[MdhBlock(iBlock).length nScanToRead],'uchar=>uchar');
  
      offset = MDH_SCANSIZE + MDH_PMUSIZE;

      % PMU Header
      MdhParam = extract_mdh_parameters(blockOfScans(offset+1:end,:), 'pmu', {'timeStamp', 'duration', 'data'});
      pmuTimeStamp = reshape(MdhParam.timeStamp,2,nScanToRead);
      pmuDuration = MdhParam.duration;

      PMU.TimeStamp.raw(:,PMU.TimeStamp.nColumn+1:PMU.TimeStamp.nColumn + nScanToRead) = pmuTimeStamp;
      PMU.TimeStamp.data(:,PMU.TimeStamp.nColumn+1:PMU.TimeStamp.nColumn + nScanToRead) = double(pmuTimeStamp) * 2.5d-3;
      PMU.TimeStamp.nColumn = PMU.TimeStamp.nColumn + nScanToRead;

      pmuDataOld = blockOfScans(offset+17:end,:);
      nPmuPoints = idivide(MdhBlock(iBlock).length - offset - uint64(16),uint64(4),'fix');
      pmuData = reshape(MdhParam.data,nPmuPoints,nScanToRead);

      iPmuPoint = 1;
      while iPmuPoint <= nPmuPoints

         pmuType = bitshift(bitand(1.0*pmuData(iPmuPoint,:),hex2dec('010F0000'),'uint32'),-16) - 256;
         if length(unique(pmuType)) ~= 1
            if ~Control.User.silent
               display_message(Control.User.noGui,'ERROR','Unexpected PMU type data','PMU Problem');
            end
            fclose(fid);
            PMUOut = [];
            return;
         else
            pmuType = unique(pmuType);
         end
         % Check for end of data magic
         if (pmuType < 1) | (pmuType > 8)
            break 
         end
         pmuPeriod = pmuData(iPmuPoint+1,:);
         nPeriod = idivide(pmuDuration(:),pmuPeriod(:),'fix');
         if (length(unique(nPeriod)) ~= 1) || (length(unique(pmuPeriod(:))) ~= 1)
            if ~Control.User.silent
               display_message(Control.User.noGui,'ERROR','Unexpected PMU period or point data','PMU Problem');
            end
            fclose(fid);
            PMUOut = [];
            return;
         else
            nPeriod = unique(nPeriod);
            pmuPeriod = unique(pmuPeriod);
         end
         
         pmuPoint = uint32(iPmuPoint);
         pmuThisData = pmuData(pmuPoint+2:pmuPoint+2+nPeriod-1,:);
         % Waveforms are the lower 12 bits
         pmuWaveform = bitand(1.0*pmuThisData,hex2dec('00000fff'),'uint32');
         % Triggers in higher bits
         pmuTriggers = bitand(1.0*pmuThisData,hex2dec('fffff000'),'uint32');

         if PMU.Waveforms(pmuType).period == 0
            PMU.Waveforms(pmuType).period = (double(pmuPeriod) * 100d-6);
         else
            if PMU.Waveforms(pmuType).period ~= (double(pmuPeriod)*100d-6)
               if ~Control.User.silent
                  display_message(Control.User.noGui,'ERROR','Unexpected PMU period','PMU Problem');
               end
               fclose(fid);
               PMUOut = [];
               return;

            end
         end
         PMU.Waveforms(pmuType).data(1:nPeriod,PMU.Waveforms(pmuType).nColumn+1:PMU.Waveforms(pmuType).nColumn + nScanToRead) = pmuWaveform;
         PMU.Waveforms(pmuType).trigger(1:nPeriod,PMU.Waveforms(pmuType).nColumn+1:PMU.Waveforms(pmuType).nColumn + nScanToRead) = pmuTriggers;
         PMU.Waveforms(pmuType).nRow(PMU.Waveforms(pmuType).nColumn+1:PMU.Waveforms(pmuType).nColumn + nScanToRead) = nPeriod;
         PMU.Waveforms(pmuType).nColumn = PMU.Waveforms(pmuType).nColumn + nScanToRead;

         % Need to jump our pointer forward
         iPmuPoint = iPmuPoint + 2 + double(nPeriod);

      end

   end % for over blocks

   % Consolidate values
   WaveformOut.typeName = '';
   WaveformOut.nPoint = 0;
   WaveformOut.time = [];
   WaveformOut.data = [];
   WaveformOut.trigger = [];
   pmuTypes = {'ECG1', 'ECG2', 'ECG3', 'ECG4', 'PULS', 'RESP', 'EXT1', 'EXT2'};
   PMUOut.TimeStamp = PMU.TimeStamp;
   PMUOut.Waveforms = repmat(WaveformOut,1,length(pmuTypes));
   for iType=1:length(pmuTypes)
      PMUOut.Waveforms(iType).typeName = pmuTypes(iType);
      nRowTotal = sum(PMU.Waveforms(iType).nRow(:));
      PMUOut.Waveforms(iType).data = uint32(0) * zeros(1,nRowTotal,'uint32');
      PMUOut.Waveforms(iType).trigger = uint32(0)*zeros(1,nRowTotal,'uint32');
      PMUOut.Waveforms(iType).time = 0 * zeros(1,nRowTotal);
      iStart = 1;
      PMUOut.Waveforms(iType).nPoint = nRowTotal;
      for iCol=1:PMU.Waveforms(iType).nColumn

         nRow = PMU.Waveforms(iType).nRow(iCol);
         timeRelative = (0:(nRow-1)) * PMU.Waveforms(iType).period;
         timeStart = PMU.TimeStamp.data(2,iCol);
         timeValues = timeStart + timeRelative;
         waveValues = PMU.Waveforms(iType).data(1:nRow,iCol);
         triggers = PMU.Waveforms(iType).trigger(1:nRow,iCol);

         PMUOut.Waveforms(iType).time(iStart:iStart+nRow-1) = timeValues(:);
         PMUOut.Waveforms(iType).data(iStart:iStart+nRow-1) = waveValues(:);
         PMUOut.Waveforms(iType).trigger(iStart:iStart+nRow-1) = triggers(:);

         iStart = iStart + nRow;

      end
   end

   clear PMU;

end

% =============================================================================
function [dimNames, nScan, n, kSpaceCentreLineNo, kSpaceCentrePartitionNo, ...
          order, dimOut, dimOutProd, channelIDUniq] = ...
         scan_mdh_for_dim_and_order(fid,Dim,Control,sizeChanInBytes,MdhBlock,nScanInMeas,nSamplesInScan)
% -----------------------------------------------------------------------------
% -----------------------------------------------------------------------------

   global MDH_SCANSIZE MDH_CHANSIZE MDH_PMUSIZE MAX_BYTES_PER_FREAD

   iScanRead = uint64(0);
   dimNames  = {'CHA', 'LIN', 'SLC', 'PAR', 'ACQ', 'ECO', 'PHS', 'REP', 'SET'};
   indexPrev.dumby = 0; % Placeholder
   for typeCell = {'PhaseCor', 'RTFeedback', 'NoiseAdj', 'Normal'}
      dimValue.(typeCell{1}) = uint16(0)*uint16(zeros(1,length(dimNames)));
      dimRank.(typeCell{1}) = 0 * zeros(1,length(dimNames));
      nScan.(typeCell{1}) = 0;
      n.SEG.(typeCell{1}) = 0;
   end
   if ~(Control.User.noGui || Control.User.silent)
      time0 = tic;
   end

   kSpaceCentreLineNo = [];
   kSpaceCentrePartitionNo = [];

   % Instead of looping through blocks, we partition them into meta-blocks
   %    of the same number of samples per scan
   iBlock = 1;
   iBlockMeta = 1;
   MdhBlockMeta = append_to_mdhblockmeta(iBlock, iBlockMeta, MdhBlock);
   for iBlock=2:length(MdhBlock)

      indexMatch = find((cell2mat({MdhBlockMeta.isSyncData}) == MdhBlock(iBlock).isSyncData) & ...
                        (cell2mat({MdhBlockMeta.isAcqEnd}) == MdhBlock(iBlock).isAcqEnd) & ...
                        (cell2mat({MdhBlockMeta.length}) == MdhBlock(iBlock).length) & ...
                        (cell2mat({MdhBlockMeta.nSamplesInScan}) == MdhBlock(iBlock).nSamplesInScan) & ...
                        (cell2mat({MdhBlockMeta.nChannelUsed}) == MdhBlock(iBlock).nChannelUsed));
      if numel(indexMatch) > 0
         MdhBlockMeta(indexMatch(1)).nBlock = MdhBlockMeta(indexMatch(1)).nBlock + 1;
         MdhBlockMeta(indexMatch(1)).iBlock(MdhBlockMeta(indexMatch(1)).nBlock) = iBlock;
         MdhBlockMeta(indexMatch(1)).nScan = MdhBlockMeta(indexMatch(1)).nScan + MdhBlock(iBlock).nScan;
      else 
         iBlockMeta = iBlockMeta + 1;
         MdhBlockMeta = append_to_mdhblockmeta(iBlock, iBlockMeta, MdhBlock, MdhBlockMeta);
      end
   end % for over MdhBlock

   for iBlockMeta = 1:length(MdhBlockMeta)

      % Skip SYNCDATA
      if MdhBlockMeta(iBlockMeta).isSyncData || MdhBlockMeta(iBlockMeta).isAcqEnd
         continue
      end

      for iBlockInner = 1:length(MdhBlockMeta(iBlockMeta).iBlock)

         iBlock = MdhBlockMeta(iBlockMeta).iBlock(iBlockInner);

         maxScanPerRead = max(idivide(MAX_BYTES_PER_FREAD,MdhBlock(iBlock).length,'fix'),uint64(1));
         nScanToRead = min(maxScanPerRead,MdhBlock(iBlock).nScan);

         fseek(fid,MdhBlock(iBlock).offsetInFile,'bof');
         nScanRemainInBlock = MdhBlock(iBlock).nScan;

         while nScanRemainInBlock > 0

            if ~(Control.User.noGui || Control.User.silent)
               dTimeElapsed = toc(time0);
               rate = dTimeElapsed / (double(iScanRead) + 1);
               timeRemain = double(nScanInMeas - iScanRead) * rate;
               textMessage = sprintf('Pass 2 of 3\nScanning dimensions from MDH\nTime remaining (s): %d', uint32(timeRemain));
               waitbar(double(iScanRead)/double(nScanInMeas), Control.GUI.hWaitBar, textMessage);
            end

            blockOfScans = fread(fid,[MdhBlock(iBlock).length nScanToRead],'uchar=>uchar');
            MdhParam = extract_mdh_parameters(blockOfScans, 'scan', {'indexLIN','indexACQ','indexSLC','indexPAR','indexECO','indexPHS','indexREP','indexSET','indexSEG'});
            MdhBitFlags = extract_mdh_bitflags(blockOfScans, {'isAcqEnd', 'isSyncData', 'isNoiseAdj', 'isPhaseCor', 'isRTFeedback', 'isRawDataCorrection', 'isReflect'});

            bHasACQEND = max(MdhBitFlags.isAcqEnd);
            if bHasACQEND
               break
            end
            channelID = uint16(0) * zeros(MdhBlock(iBlock).nChannelUsed,nScanToRead,'uint16');
            for iCha = 1:MdhBlock(iBlock).nChannelUsed
               offset  = uint64(MDH_SCANSIZE) + uint64(iCha-1) * sizeChanInBytes;
               channelID(iCha,:) = reshape(extract_single_mdh_parameter(blockOfScans(offset+1:offset+sizeChanInBytes,:),'chan','channelID'),[1,nScanToRead]);
            end
            channelIDUniq = unique(channelID(:));
            channelIDUniqNum = length(channelIDUniq);
            if (channelIDUniqNum ~= MdhBlock(iBlock).nChannelUsed) || (sum(channelIDUniq == channelID(:,1)) ~= MdhBlock(iBlock).nChannelUsed)
               if ~Control.User.silent
                  display_message(Control.User.noGui,'ERROR',['Mismatch between nChannel and coils found: ' int2str(channelIDUniqNum) ' ' int2str(MdhBlock(iBlock).nChannelUsed)],'Dimension Conflict');
               end
               fclose(fid);
               dimNames = [];
               nScan = [];
               n = [];
               kSpaceCentreLineNo = [];
               kSpaceCentrePartitionNo = [];
               order = [];
               dimOut = [];
               dimOutProd = [];
               channelIDUniq = [];
               return;
            end
            maskNormal = MdhBitFlags.isNoiseAdj == 0;
            for typeCell = {'NoiseAdj', 'PhaseCor', 'RTFeedback','Normal'}
               typeStr = typeCell{1};
               bitField = strcat('is',typeStr);
               readField = strcat('read',typeStr);
               if strcmp(typeStr,'Normal')
                  index.(typeStr) = find(maskNormal);
               else
                  index.(typeStr) = find(MdhBitFlags.(bitField));
                  maskNormal = maskNormal & (MdhBitFlags.(bitField) == 0);
               end
               if Control.User.(readField) && (numel(index.(typeStr)) > 0)
                  if strcmp(typeStr,'Normal')
                     if nScan.Normal == 0
                        kSpaceCentreLineNo = extract_single_mdh_parameter(blockOfScans(:,index.Normal(1)),'scan','kSpaceCentreLineNo');
                        kSpaceCentrePartitionNo = extract_single_mdh_parameter(blockOfScans(:,index.Normal(1)),'scan','kSpaceCentrePartitionNo');
                     end
                  end
                  nScan.(typeStr) = nScan.(typeStr) + numel(index.(typeStr));
                  n.SEG.(typeStr) = max(n.SEG.(typeStr), max(MdhParam.indexSEG(index.(typeStr)(:))));
                  dimValue.(typeStr) = max_from_mdh_param(dimValue.(typeStr), MdhBlock(iBlock).nChannelUsed, dimNames, MdhParam, index.(typeStr)(:));
                  if numel(index.(typeStr)) > 1
                     dimRank.(typeStr) = rank_from_mdh_param(dimRank.(typeStr), dimNames, MdhParam, MdhBlock(iBlock).nChannelUsed, index.(typeStr)(:));
                  end
                  if numel(index.(typeStr)) == 1
                     dimRank.(typeStr)(1) = dimRank.(typeStr)(1) + MdhBlock(iBlock).nChannelUsed;
                  end
                  if (iBlockInner > 1) && is_field_in_struct(indexPrev, typeStr);
                     dimRank.(typeStr) = rank_from_index_prev(dimRank.(typeStr), ...
                                                                dimNames, ...
                                                                MdhParam, ...
                                                                MdhBlock(iBlock).nChannelUsed, ...
                                                                index.(typeStr)(1), ...
                                                                indexPrev.(typeStr));
                  end
                  if MdhBlock(iBlock).nChannelUsed < 2
                     dimRank.(typeStr)(1)=0;
                  end
               end
               if numel(index.(typeStr)) > 0
                  indexPrev.(typeStr) = set_index_prev(dimNames, MdhParam, index.(typeStr)(end));
               else
                  if is_field_in_struct(indexPrev,typeStr)
                     indexPrev = rmfield(indexPrev,typeStr);
                  end
               end
            end

            iScanRead = iScanRead + nScanToRead;
            nScanRemainInBlock = MdhBlock(iBlock).nScan - iScanRead;
            nScanToRead = min(nScanToRead,nScanRemainInBlock);

         end % while over subsets of blocked lines

      end % for over blocks of non-syncdata

   end % over blocks of the same size

   if ~(Control.User.noGui || Control.User.silent)
      waitbar(0.5, Control.GUI.hWaitBar, 'MDH Scan complete');
   end
   clear blockOfScans;
   clear channelID MdhParam MdhBitFlags;
   clear index;
   clear iScanRead nScanToRead;

   for typeCell = {'PhaseCor', 'NoiseAdj', 'RTFeedback', 'Normal'}
      typeStr = typeCell{1};
      dimValue.(typeStr) = dimValue.(typeStr) .* uint16(dimRank.(typeStr) > 0) + 1;
      for dimCell = {{'CHA' 1 false ''    ''     false ''      ''                   ''      ''    } ...
                     {'LIN' 2 false ''    ''     true  'Recon' 'nFourierLines'      'Recon' 'nLin'} ...
                     {'SLC' 3 false ''    ''     true  ''      ''                   'Recon' 'nSlc'} ...
                     {'PAR' 4 false ''    ''     true  'Recon' 'nFourierPartitions' 'Recon' 'nPar'} ...
                     {'ACQ' 5 true  'Raw' 'nAcq' false ''      ''                   ''      ''    } ...
                     {'ECO' 6 true  'Raw' 'nEco' false ''      ''                   ''      ''    } ...
                     {'PHS' 7 true  'Raw' 'nPhs' false ''      ''                   ''      ''    } ...
                     {'REP' 8 true  'Raw' 'nRep' false ''      ''                   ''      ''    } ...
                     {'SET' 9 true  'Raw' 'nSet' false ''      ''                   ''      ''    }}
         % Initialize with dimValue
         n.(dimCell{1}{1}).(typeStr) = dimValue.(typeStr)(dimCell{1}{2});
         % Update with max from Dim structure
         if dimCell{1}{3} & (n.(dimCell{1}{1}).(typeStr) > 1)
            n.(dimCell{1}{1}).(typeStr) = max(n.(dimCell{1}{1}).(typeStr), Dim.(dimCell{1}{4}).(dimCell{1}{5}));
         end
         % Special handling in FillToFullFourier cases
         if dimCell{1}{6}
            if ~Control.User.noFillToFullFourier
               if length(dimCell{1}{7}) > 0
                  n.(dimCell{1}{1}).(typeStr) = max(n.(dimCell{1}{1}).(typeStr), Dim.(dimCell{1}{7}).(dimCell{1}{8}));
               end
               if length(dimCell{1}{9}) > 0
                  n.(dimCell{1}{1}).(typeStr) = max(n.(dimCell{1}{1}).(typeStr), Dim.(dimCell{1}{9}).(dimCell{1}{10}));
               end
            end
         end
      end
      n.SEG.(typeStr) = n.SEG.(typeStr) + 1;
      dimOut.(typeStr) = uint64([n.CHA.(typeStr) n.LIN.(typeStr) n.SLC.(typeStr) n.PAR.(typeStr) n.ACQ.(typeStr) n.ECO.(typeStr) n.PHS.(typeStr) n.REP.(typeStr) n.SET.(typeStr)]);
      [~, order.(typeStr)] = sort(dimRank.(typeStr), 2, 'descend');
      dimOut.(typeStr) = dimOut.(typeStr)(order.(typeStr));
      dimOutProd.(typeStr) = [uint64(1) compute_cumulative_product(dimOut.(typeStr)(1:end-1))];

   end
   if Control.User.noFillToFullFourier == true 
      n.COL = nSamplesInScan;
   else
      n.COL = max(nSamplesInScan, uint64(Dim.Recon.readoutFTLength));
   end

   clear dimValue dimRank;

end

% =============================================================================
% -----------------------------------------------------------------------------
% -----------------------------------------------------------------------------

