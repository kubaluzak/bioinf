import xml.etree.ElementTree as ET
import os
import time
from collections import deque, Counter
import heapq

def process_file(filepath, beam_width=2):
    root = ET.parse(filepath).getroot()
    start_seq = root.attrib['start']
    target_length = int(root.attrib['length'])
    probe_elements = root.findall('probe')
    k = len(probe_elements[0].attrib['pattern']) - 1

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
        return ''.join({'A': 'W', 'T': 'W', 'C': 'S', 'G': 'S'}[n] for n in seq)

    def encode_ry(seq):
        return ''.join({'A': 'R', 'G': 'R', 'C': 'Y', 'T': 'Y'}[n] for n in seq)

    # Dodatkowa heurystyka: ranking liter wg częstotliwości występowania sond
    freq_counter = Counter()
    for p in probes1 + probes2:
        freq_counter[p[-1]] += 1

    start_suffix_ws = encode_ws(start_seq)[-k:]
    start_suffix_ry = encode_ry(start_seq)[-k:]

    visited = set()
    heap = []

    # Heap -> (priorytet, długość seq, seq, s_ws, s_ry)
    heapq.heappush(heap, (-0, len(start_seq), start_seq, start_suffix_ws, start_suffix_ry))
    visited.add((start_suffix_ws, start_suffix_ry))

    result = None

    start_time = time.perf_counter()

    while heap:
        new_heap = []
        for _ in range(min(len(heap), beam_width)):
            _, _, seq, s_ws, s_ry = heapq.heappop(heap)

            if len(seq) >= target_length:
                result = seq
                end_time = time.perf_counter() - start_time
                return k, end_time, len(result), result

            cands_1 = probes_map_1.get(s_ws, set())
            cands_2 = probes_map_2.get(s_ry, set())

            for N in ['A', 'C', 'G', 'T']:
                if N in cands_1 and N in cands_2:
                    new_s_ws = s_ws[1:] + encode_ws(N)
                    new_s_ry = s_ry[1:] + encode_ry(N)
                    if (new_s_ws, new_s_ry) not in visited:
                        visited.add((new_s_ws, new_s_ry))
                        # Priorytet = częstotliwość litery + długość sekwencji (heurystyka)
                        priority = freq_counter[N] + len(seq)
                        heapq.heappush(new_heap, (-priority, len(seq) + 1, seq + N, new_s_ws, new_s_ry))

        heap = new_heap

    end_time = time.perf_counter() - start_time
    seq_len = len(result) if result else 0
    return k, end_time, seq_len, result


def process_folder(folder_path):
    files = [f for f in os.listdir(folder_path) if f.lower().endswith('.xml')]
    files = sorted(files)

    results = []
    for filename in files:
        filepath = os.path.join(folder_path, filename)
        k, elapsed, seq_len, result = process_file(filepath, beam_width=2)
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
