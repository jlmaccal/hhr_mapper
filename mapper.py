def find_line_ranges(lines):
    start = None
    ranges = []
    for i, line in enumerate(lines):
        if line.startswith('No'):
            cols = line.split()
            # check to see we have 'No' followed by a number and nothing else
            if not len(cols) == 2:
                continue
            try:
                dummy = int(cols[1])
            except:
                continue

            # we've found the start of a line
            # if we were previously dealing with a line, we need to store it
            if not start is None:
                ranges.append((start, i))

            # now we set start to the current line and search again
            start = i

        if line.startswith('Done!'):
            # this is the last line of the file
            # we'll store everything up to the 'Done!' as the last line
            ranges.append((start, i))

    return ranges


def parse_match_number(line):
    assert line.startswith('No')
    cols = line.split()
    return int(cols[1])


def parse_description(line):
    assert line.startswith('>')
    cols = line.split()
    pdb_string = cols[0][1:]
    pdb, chain = pdb_string.split('_')
    return pdb, chain, ' '.join(cols[1:])


def get_sequence_lines(lines, ident, skip_ss=True):
    sequence_lines = []
    for line in lines:
        if line.startswith('{} '.format(ident)):
            if not ('ss_pred' in line or
                    'Consensus' in line or
                    'ss_dssp' in line or
                    'ss_pred' in line) or not skip_ss:
                sequence_lines.append(line)
    return sequence_lines


def get_query_sequence_lines(lines):
    return get_sequence_lines(lines, 'Q')


def get_template_sequence_lines(lines):
    return get_sequence_lines(lines, 'T')


def get_query_ss_pred_lines(lines):
    return get_sequence_lines(lines, 'Q ss_pred', skip_ss=False)


def get_template_ss_pred_lines(lines):
    return get_sequence_lines(lines, 'T ss_pred', skip_ss=False)


def get_template_dssp_lines(lines):
    return get_sequence_lines(lines, 'T ss_dssp', skip_ss=False)


def extract_residue_range_from_sequence_line(line):
    cols = line.split()
    start = int(cols[2])
    end = int(cols[4])
    return start, end


def extract_sequence_from_sequence_line(line):
    cols = line.split()
    return cols[3]


def extract_ss_from_ss_line(line):
    cols = line.split()
    return cols[2]


def join_ss_lines(ss_lines):
    ss = ''.join([extract_ss_from_ss_line(ss) for ss in ss_lines])
    return ss


def get_query_ss(lines):
    ss_lines = get_query_ss_pred_lines(lines)
    return join_ss_lines(ss_lines)


def get_template_ss(lines):
    ss_lines = get_template_ss_pred_lines(lines)
    return join_ss_lines(ss_lines)


def get_template_dssp(lines):
    ss_lines = get_template_dssp_lines(lines)
    return join_ss_lines(ss_lines)


def join_seq(seq_lines):
    seq = ''.join([extract_sequence_from_sequence_line(s) for s in seq_lines])
    numbers = [extract_residue_range_from_sequence_line(s) for s in seq_lines]
    start = numbers[0][0]
    end = numbers[-1][1]
    n = len(numbers)
    for i in range(n - 1):
        assert numbers[i][1] + 1 == numbers[i + 1][0]
    return seq, start, end


def get_query_sequence(lines):
    seq_lines = get_query_sequence_lines(lines)
    return join_seq(seq_lines)


def get_template_sequence(lines):
    seq_lines = get_template_sequence_lines(lines)
    return join_seq(seq_lines)
