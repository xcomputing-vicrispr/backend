def getMMsequence_flexible(original, bowtie_details, strand="+"):


    entries = [x.strip() for x in bowtie_details.split(";") if x.strip()]

    if not original or not bowtie_details:
        return "---------"
        
    original = original.upper()
    seq_list_original = list(original.upper())
    final_res = []
    final_pos = []
    
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

    try:
        for entry in entries: 
            pos = []
            if ":" not in entry:
                return ""
            
            if strand == "-":
                working_seq = "".join(complement_map.get(base, base) for base in entry)
            else:
                working_seq = entry

            rules = working_seq.split(",")
            
            for rule in rules:
                if ":" not in rule or ">" not in rule:
                    continue
                    
                pos_part, change_part = rule.split(":")
                idx = int(pos_part)
                new_char, old_char = change_part.split(">")

                if idx < 0 or idx >= len(original) or original[idx] != old_char:
                    return "-------"
                
                seq_list_original[idx] = new_char
                pos.append(idx)
            
            res = "".join(seq_list_original)
            final_res.append(res)
            final_pos.append(pos)

        
    except (ValueError, IndexError, KeyError) as e:
        print(f"Mapping error: {e}")
        return "---------"

    return final_res, final_pos

x = getMMsequence_flexible("AATAAAAAAATTTTATCTTGTGG", "CP019943.1:46490,,15:T>A,17:T>A,18:C>A;", "-")

print(x)

