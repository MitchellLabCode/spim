function combine_xml_files(xml_files, output_file)
    % Parse multiple XML files and merge ViewInterestPointsFile nodes while preserving structure
    % xml_files: Cell array of XML file paths
    % output_file: Output XML file path
    
    % Use the first XML file as a template
    base_doc = xmlread(xml_files{1});
    root_elem = base_doc.getDocumentElement();
    
    % Locate ViewInterestPoints section
    vp_nodes = base_doc.getElementsByTagName('ViewInterestPoints');
    if vp_nodes.getLength == 0
        error('No <ViewInterestPoints> section found in template XML.');
    end
    vp_elem = vp_nodes.item(0);
    
    % Clear existing ViewInterestPointsFile entries
    existing_files = vp_elem.getElementsByTagName('ViewInterestPointsFile');
    while existing_files.getLength > 0
        vp_elem.removeChild(existing_files.item(0));
    end
    
    % Collect and merge data from all XML files
    data = struct();
    keys = {};
    
    for i = 1:length(xml_files)
        xml_doc = xmlread(xml_files{i});
        nodes = xml_doc.getElementsByTagName('ViewInterestPointsFile');
        
        for j = 0:nodes.getLength-1
            node = nodes.item(j);
            timepoint = str2double(node.getAttribute('timepoint'));
            setup = str2double(node.getAttribute('setup'));
            params = strtrim(char(node.getAttribute('params'))); % Trim whitespace
            content = strtrim(char(node.getTextContent())); % Trim whitespace
            
            key = sprintf('T%d_S%d', timepoint, setup); % Ensure valid struct field name
            
            if isfield(data, key)
                % Check for conflicting params
                if ~strcmp(data.(key).params, params)
                    error('Conflicting params for timepoint=%d, setup=%d', timepoint, setup);
                end
            else
                keys{end+1} = key; %#ok<AGROW>
            end
            
            data.(key) = struct('timepoint', timepoint, 'setup', setup, 'params', params, 'content', content);
        end
    end
    
    % Extract values and sort by timepoint and setup
    num_entries = numel(keys);
    numeric_data = zeros(num_entries, 2);
    
    for i = 1:num_entries
        numeric_data(i, :) = [data.(keys{i}).timepoint, data.(keys{i}).setup];
    end
    
    [~, sort_idx] = sortrows(numeric_data, [1, 2]);
    sorted_keys = keys(sort_idx);
    
    % Insert sorted ViewInterestPointsFile elements back into XML
    for i = 1:num_entries
        entry = data.(sorted_keys{i});
        elem = base_doc.createElement('ViewInterestPointsFile');
        elem.setAttribute('timepoint', num2str(entry.timepoint));
        elem.setAttribute('setup', num2str(entry.setup));
        elem.setAttribute('label', 'beads');
        elem.setAttribute('params', entry.params);
        elem.appendChild(base_doc.createTextNode(entry.content));
        vp_elem.appendChild(elem);
    end
    
    % Save to output file with pretty formatting and remove excess whitespace
    save_xml_manually(base_doc, output_file);
    fprintf('Merged XML saved to %s\n', output_file);
end
% 
% function save_xml_manually(doc, output_file)
%     % Save XML using pretty formatting by setting a document serializer
%     serializer = javax.xml.transform.stream.StreamResult(output_file);
%     transformer = javax.xml.transform.TransformerFactory.newInstance().newTransformer();
% 
%     % Enable pretty printing
%     transformer.setOutputProperty(javax.xml.transform.OutputKeys.INDENT, "yes");
%     transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "2");
% 
%     % Write the formatted XML to the file
%     transformer.transform(javax.xml.transform.dom.DOMSource(doc), serializer);
% 
%     % Now, manually remove excess whitespace between tags
%     % Read the file, remove extra whitespace, and overwrite the file
%     xml_str = fileread(output_file);
%     xml_str = regexprep(xml_str, '>\\s+<', '><'); % Remove excess spaces
%     fid = fopen(output_file, 'w');
%     fprintf(fid, '%s', xml_str);
%     fclose(fid);
% end


function save_xml_manually(doc, output_file)
    % Save XML using pretty formatting by setting a document serializer
    serializer = javax.xml.transform.stream.StreamResult(output_file);
    transformer = javax.xml.transform.TransformerFactory.newInstance().newTransformer();
    
    % Enable pretty printing
    transformer.setOutputProperty(javax.xml.transform.OutputKeys.INDENT, "yes");
    transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "2");
    
    % Write the formatted XML to the file
    transformer.transform(javax.xml.transform.dom.DOMSource(doc), serializer);
    
    % Now, manually remove excess whitespace between tags and extra blank lines
    xml_str = fileread(output_file);
    
    % Remove excess spaces between tags
    xml_str = regexprep(xml_str, '>\\s+<', '><');
    
    % Remove extra blank lines (collapse multiple newlines into a single one)
    xml_str = regexprep(xml_str, '\n\s*\n', '\n');
    
    % Write the cleaned XML back to the file
    fid = fopen(output_file, 'w');
    fprintf(fid, '%s', xml_str);
    fclose(fid);
end
