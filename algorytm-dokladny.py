import xml.etree.ElementTree as ET
from collections import deque
import os
import time

def process_file(filepath):
    root = ET.parse(filepath).getroot()
    start_seq = root.attrib['start']
    target_length = int(root.attrib['length'])
    probe_elements = root.findall('probe')
    pattern1 = probe_elements[0].attrib['pattern']
    k = len(pattern1) - 1

    probes1 = [cell.text for cell in probe_elements[0].findall('cell')]
    probes2 = [cell.text for cell in probe_elements[1].findall('cell')]

    def build_prefix_map(probes):
        m = {}
        for p in probes:
            prefix = p[:-1]
            N = p[-1]
            if prefix not in m:
                m[prefix] = set()
            m[prefix].add(N)
        return m

    probes_map_1 = build_prefix_map(probes1)
    probes_map_2 = build_prefix_map(probes2)

    def encode_ws(seq):
        map_ws = {'A':'W','T':'W','C':'S','G':'S'}
        return ''.join(map_ws[n] for n in seq)

    def encode_ry(seq):
        map_ry = {'A':'R','G':'R','C':'Y','T':'Y'}
        return ''.join(map_ry[n] for n in seq)

    start_suffix_ws = encode_ws(start_seq)[-k:]
    start_suffix_ry = encode_ry(start_seq)[-k:]

    queue = deque()
    visited = set()

    queue.append((start_suffix_ws, start_suffix_ry, start_seq))
    visited.add((start_suffix_ws, start_suffix_ry))

    result = None

    start_time = time.perf_counter()

    while queue:
        s_ws, s_ry, seq = queue.popleft()

        if len(seq) >= target_length:
            result = seq
            break

        cands_1 = probes_map_1.get(s_ws, set())
        cands_2 = probes_map_2.get(s_ry, set())

        for N in ['A', 'C', 'T', 'G']:
            if N in cands_1 and N in cands_2:
                N_ws = encode_ws(N)
                N_ry = encode_ry(N)
                new_s_ws = s_ws[1:] + N_ws
                new_s_ry = s_ry[1:] + N_ry
                if (new_s_ws, new_s_ry) not in visited:
                    visited.add((new_s_ws, new_s_ry))
                    queue.append((new_s_ws, new_s_ry, seq + N))

    end_time = time.perf_counter() - start_time

    seq_len = len(result) if result else 0

    return k, end_time, seq_len, result

def process_folder(folder_path):
    files = [f for f in os.listdir(folder_path) if f.lower().endswith('.xml')]
    files = sorted(files)

    results = []

    for filename in files:
        filepath = os.path.join(folder_path, filename)
        k, elapsed, seq_len, result = process_file(filepath)
        results.append((filename, k, elapsed, seq_len, result))

    results.sort(key=lambda x: x[3])

    print(f"{'Długość sekwencji':<18} {'Czas [s]':<10}  {'Sekwencja'}")
    print("-" * 100)

    for filename, k, elapsed, seq_len, result in results:
        status = result if result else "BRAK ROZWIĄZANIA"
        print(f"{seq_len:<18} {elapsed:<10.6f}  {status}")

if __name__ == "__main__":
    folder = "bio_files"  # podaj ścieżkę do folderu z plikami XML
    process_folder(folder)
