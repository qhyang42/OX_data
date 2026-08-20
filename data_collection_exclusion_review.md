# Data collection notes: potential exclusions and QC checks

## Scope

This is a notes-only review of the handwritten PDFs in `data_collection_notes/` for subjects 2-6. It records events that may justify excluding a trial/run or replacing an aborted run. It is not a final exclusion decision. Confirm each item against behavioral output, acquisition logs, respiratory/lab-chart traces, and motion estimates before changing the analysis dataset.

Interpretation used below:

- **Exclude/replace**: the note explicitly says a run was aborted, incomplete, redone, or had no/missed data.
- **Review**: a discrete motion, cough, weak-odor, or behavioral event was noted; apply the project's quantitative QC rule before exclusion.
- **Session QC**: the note identifies a possible problem but does not provide enough trial-level detail.

## Subject 2

### Trial/run candidates

| Session | Run | Trial | Note-derived issue | Preliminary handling | Source note |
|---|---:|---:|---|---|---|
| 7 | 3 | 9 | "Possible cough" | Review trial against motion and respiratory trace. | [session 7](data_collection_notes/subj2/240723_fMRI_OX_NWU_LS_note_session7.pdf) |
| 9 | 1 | 6 | No odor. | Exclude the affected odor trial unless delivery logs contradict the note. | [sessions 8 and 9, second note sheet](data_collection_notes/subj2/240723_fMRI_OX_NWU_LS_note_session8and9.pdf) |

### Session-level QC

| Session | Run/trial | Note-derived issue | Follow-up | Source note |
|---|---|---|---|---|
| 2 | Run 8 / unspecified runs | Participant tired after run 8; some runs had lower-amplitude breathing signal. | Check respiratory trace quality and whether fatigue affected later responses. No trial exclusion can be assigned from the note alone. | [session 2](data_collection_notes/subj2/240723_fMRI_OX_NWU_LS_note_session2.pdf) |
| 3 | All runs / unspecified | Vacuum connector was changed; odors were described as weaker afterward. | Compare odor-delivery/lab-chart traces across the session. | [session 3](data_collection_notes/subj2/240723_fMRI_OX_NWU_LS_note_session3.pdf) |
| 4 | Unspecified trials | Participant reported many no-odor trials; "MFC2 not working - check trials." | Identify affected trials from MFC/DAQ and lab-chart logs, then exclude confirmed delivery failures. | [session 4](data_collection_notes/subj2/240723_fMRI_OX_NWU_LS_note_session4.pdf) |
| 6 | Unspecified | Task marked not completed; participant was unable to concentrate. | Establish the last valid completed run/trial from behavioral logs before analysis. | [session 6](data_collection_notes/subj2/240723_fMRI_OX_NWU_LS_note_session6.pdf) |
| 7 | Imaging scan, not an OX trial | First whole-brain scan had too much motion; note says to use the second scan. | Use the second whole-brain scan and do not treat this as an OX trial exclusion. | [session 7](data_collection_notes/subj2/240723_fMRI_OX_NWU_LS_note_session7.pdf) |
| 10 | Positioning/check scan | Screen error said some documents were rejected/could not be used for positioning; note says to check the scan. | Verify the selected positioning/anatomical scan and alignment. | [session 10](data_collection_notes/subj2/240723_fMRI_OX_NWU_LS_note_session10.pdf) |

No explicit exclusion-worthy event was noted for sessions 1, 5, 8, or 11.

## Subject 3

### Trial/run candidates

| Session | Run | Trial | Note-derived issue | Preliminary handling | Source note |
|---|---:|---:|---|---|---|
| 3 | 2 | 2 or 3 | Big motion spike about 45 seconds into the run; the note is uncertain between trial 2 and trial 3. | Use timestamps/motion trace to resolve the exact trial, then apply motion QC. | [session 3](data_collection_notes/subj3/240814_fMRI_OX_NWU_JN_note_session3.pdf) |
| 6 | 5 | 2 | Large motion spike. | Review/exclude under the motion criterion. | [session 6](data_collection_notes/subj3/240814_fMRI_OX_NWU_JN_note_session6.pdf) |
| 6 | 6 | 9 | Motion attributed to an itch. | Review/exclude under the motion criterion. | [session 6](data_collection_notes/subj3/240814_fMRI_OX_NWU_JN_note_session6.pdf) |
| 6 | 7-8 | All | Task stopped because of a computer problem; the note lists only runs 1-6 as collected. | Treat runs 7-8 as not collected, not as valid runs. | [session 6](data_collection_notes/subj3/240814_fMRI_OX_NWU_JN_note_session6.pdf) |
| 9 | 7 | Unspecified | Session stopped during run 7 because of ear itch. | Treat run 7 as incomplete unless logs show it finished before the stop. | [session 9](data_collection_notes/subj3/240814_fMRI_OX_NWU_JN_note_session9.pdf) |
| 9 | 8 | All | Run 8 was planned but not collected after the run-7 stop. | Treat as not collected. | [session 9](data_collection_notes/subj3/240814_fMRI_OX_NWU_JN_note_session9.pdf) |

### Session-level QC

| Session | Run/trial | Note-derived issue | Follow-up | Source note |
|---|---|---|---|---|
| 1 | Before task; boundary between runs 5-6 | Task did not start initially and was fixed by restarting the laptop; audio-device index changed between runs 5 and 6. | Verify no aborted output was retained and confirm cue/audio continuity beginning with run 6. | [session 1](data_collection_notes/subj3/240814_fMRI_OX_NWU_JN_note_session1.pdf) |
| 2 | Before run 1 | DAQ malfunction before run 1. | Verify run 1 began cleanly and contains the expected DAQ/lab-chart samples. | [session 2](data_collection_notes/subj3/240814_fMRI_OX_NWU_JN_note_session2.pdf) |
| 4 | Unspecified odor/trial in set 1 | One odor was described as very weak. | Cross-reference odor identity and delivery trace; exclude only a confirmed delivery failure. | [session 4](data_collection_notes/subj3/240814_fMRI_OX_NWU_JN_note_session4.pdf) |
| 5 | Beginning of an unspecified run/session segment | Lab-chart data showed at least three pulses at the beginning and "look[ed] weird." | Inspect trigger/pulse timing and align the affected run before trial-level decisions. | [session 5](data_collection_notes/subj3/240814_fMRI_OX_NWU_JN_note_session5.pdf) |

No explicit exclusion-worthy event was noted for sessions 7, 8, or 10.

## Subject 4

### Trial/run candidates

| Session | Run | Trial | Note-derived issue | Preliminary handling | Source note |
|---|---:|---:|---|---|---|
| 2 | 7 | Unspecified | Stopped at run 7 because of potential fatigue. | Treat run 7 as incomplete unless completion is confirmed. | [session 2](data_collection_notes/subj4/240816_fMRI_OX_NWU_RR_notes_session2.pdf) |
| 2 | 8 | All | Section/run 8 was empty after the fatigue stop. | Treat as not collected. | [session 2](data_collection_notes/subj4/240816_fMRI_OX_NWU_RR_notes_session2.pdf) |
| 3 | 6 | All/unspecified | Participant came out at run 6 for a neck adjustment and did not continue. | Treat run 6 as incomplete/not collected; runs 1-5 appear to be the completed set. | [session 3](data_collection_notes/subj4/240816_fMRI_OX_NWU_RR_notes_session3.pdf) |
| 5, aborted attempt (2024-10-23) | 1 and any additional partial data | All | Participant came out after run 1 because of neck discomfort. A later completed session-5 note exists (2024-10-29). | Exclude the aborted attempt and use the later completed acquisition after file identity is verified. | [session 5 stopped](data_collection_notes/subj4/240816_fMRI_OX_NWU_RR_notes_session5_stopped.pdf) |
| 5, completed attempt (2024-10-29) | 1 | 10 | Motion spike. | Review/exclude under the motion criterion. | [session 5](data_collection_notes/subj4/240816_fMRI_OX_NWU_RR_notes_session5.pdf) |
| 5, completed attempt (2024-10-29) | 4 | 4 | Motion spike. | Review/exclude under the motion criterion. | [session 5](data_collection_notes/subj4/240816_fMRI_OX_NWU_RR_notes_session5.pdf) |
| 5, completed attempt (2024-10-29) | 5 | 1, 4, 8 | Motion spikes. | Review/exclude affected trials under the motion criterion. | [session 5](data_collection_notes/subj4/240816_fMRI_OX_NWU_RR_notes_session5.pdf) |
| 6 | 1 | 5 | Motion spike. | Review/exclude under the motion criterion. | [session 6](data_collection_notes/subj4/240816_fMRI_OX_NWU_RR_notes_session6.pdf) |
| 8 | 4 | 3 | Motion. | Review/exclude under the motion criterion. | [session 8](data_collection_notes/subj4/240816_fMRI_OX_NWU_RR_notes_session8.pdf) |
| 11, aborted attempt | 1-5 | All | Stopped because of neck/shoulder discomfort; note says to redo the whole session next time. | Exclude the entire aborted attempt and use the later repeat, once matched. | [second note sheet inside the session-10 PDF](data_collection_notes/subj4/240816_fMRI_OX_NWU_RR_notes_session10.pdf) |

### Session/file-identity QC

| Session | Issue | Follow-up | Source note |
|---|---|---|---|
| 7 | Initially mislabeled as session 8; checklist and result name were changed to session 7. | Verify behavioral, imaging, and lab-chart filenames resolve to session 7. | [session 7](data_collection_notes/subj4/240816_fMRI_OX_NWU_RR_notes_session7.pdf) |
| 11 | The session-11 aborted note is the second page of the file named as session 10. | Do not assign this page to session 10 during provenance tracking. | [session-10 PDF containing session 11](data_collection_notes/subj4/240816_fMRI_OX_NWU_RR_notes_session10.pdf) |

No explicit exclusion-worthy event was noted for sessions 1, 4, 9, 10, or 12-15.

## Subject 5

### Trial/run candidates

| Session | Run | Trial | Note-derived issue | Preliminary handling | Source note |
|---|---:|---:|---|---|---|
| 1 | 1 | 1 | No odor. | Exclude the affected odor trial unless delivery logs contradict the note. | [session 1](data_collection_notes/subj5/241018_fMRI_OX_NWU_BN_notes_session1.pdf) |
| 5 | 1 | Unspecified | Run stopped because of a script error; subsequent files/runs are numbered 2-6. | Exclude the aborted run 1 and verify run-number mapping for the retained files. | [session 5](data_collection_notes/subj5/241018_fMRI_OX_NWU_BN_notes_session5.pdf) |
| 6 | 4 | 1 | No odor. | Exclude the affected odor trial unless delivery logs contradict the note. | [session 6](data_collection_notes/subj5/241018_fMRI_OX_NWU_BN_notes_session6.pdf) |
| 9 | 5 | All/unspecified | Participant stopped with heart racing and was unable to focus; run 5 was not completed and was to be redone. | Exclude the incomplete run 5 and use the repeat if available. | [session 9](data_collection_notes/subj5/241018_fMRI_OX_NWU_BN_notes_session9.pdf) |
| 10, aborted attempt | 1-2 | All | Stopped at run 2 when heart racing recurred; note asks whether to redo the whole session. | Keep the whole attempt out of the primary dataset pending identification of the repeat. | [session 10 stopped](data_collection_notes/subj5/241018_fMRI_OX_NWU_BN_notes_session10_stopped.pdf) |

### Session/file-identity QC

| Scanner session | Script/result label | Issue | Follow-up | Source note |
|---:|---|---|---|---|
| 13 | Results line says session 12; handwritten note says "This is 13 for scanner." | Scanner-session and result-session labels differ. | Verify mapping before merging files. | [session 13](data_collection_notes/subj5/241018_fMRI_OX_NWU_BN_notes_session13.pdf) |
| 16 | Session 15 on script/results | Note explicitly maps scanner session 16 to script session 15; runs 1-3 used the black odor box and runs 4-5 the orange odor box after a two-minute rest. Motion was noted as good/low. | Preserve the mapping; no exclusion is indicated by the note. | [session 16](data_collection_notes/subj5/241018_fMRI_OX_NWU_BN_notes_session16.pdf) |
| 17 | Session 16 on script/results | Note explicitly maps scanner session 17 to script session 16. | Preserve the mapping; no exclusion is indicated by the note. | [session 17](data_collection_notes/subj5/241018_fMRI_OX_NWU_BN_notes_session17.pdf) |

No explicit exclusion-worthy event was noted for sessions 2-4, 7-8, 11-12, or 14-17 apart from the file/session mapping checks above.

## Subject 6

### Trial/run candidates

| Session | Run | Trial | Note-derived issue | Preliminary handling | Source note |
|---|---:|---:|---|---|---|
| 3 | 1 | 5 | Weak odor. | Review odor-delivery trace and exclude if delivery failure is confirmed. | [session 3](data_collection_notes/subj6/250118_fMRI_OX_NWU_VS_notes_session3.pdf) |
| 3 | 4 | 4 | Weak odor. | Review odor-delivery trace and exclude if delivery failure is confirmed. | [session 3](data_collection_notes/subj6/250118_fMRI_OX_NWU_VS_notes_session3.pdf) |
| 4 | 6 | Trial containing odor 7 (trial number not written) | Odor 7 was "reasonably strong" but was misrated. | Use the run schedule to identify the trial, then exclude/flag the behavioral response according to the response-QC rule. | [session 4](data_collection_notes/subj6/250118_fMRI_OX_NWU_VS_notes_session4.pdf) |
| 5 | 7 | 5 | Participant forgot about the context when rating. | Exclude/flag this behavioral trial. | [session 5](data_collection_notes/subj6/250118_fMRI_OX_NWU_VS_notes_session5.pdf) |
| 6 | Old run 6 | 3 | Motion spike at response. The run was subsequently redone. | Exclude the old run 6; do not also treat its individual trials as part of the retained dataset. | [session 6](data_collection_notes/subj6/250118_fMRI_OX_NWU_VS_notes_session6and7.pdf) |
| 6 | Old run 6 | 5 | Motion. The run was subsequently redone. | Exclude the old run 6. | [session 6](data_collection_notes/subj6/250118_fMRI_OX_NWU_VS_notes_session6and7.pdf) |
| 6 | Replacement mapping | All | Note says "use run 1 2 3 4 5 8," "redo run 6," and "remove section 7 (old run 6)." | Retain runs 1-5 and run 8 as the replacement for run 6; remove the old run-6/section-7 data after verifying filenames. | [session 6](data_collection_notes/subj6/250118_fMRI_OX_NWU_VS_notes_session6and7.pdf) |
| 8 | 1 | 2 | Motion spike. | Review/exclude under the motion criterion. | [session 8](data_collection_notes/subj6/250118_fMRI_OX_NWU_VS_notes_session8.pdf) |
| 8 | 2 | 3 | Potentially missed intensity response. | Confirm in behavioral output; exclude/flag the missing-response component or trial per analysis policy. | [session 8](data_collection_notes/subj6/250118_fMRI_OX_NWU_VS_notes_session8.pdf) |
| 8 | 3 | 8-9 | Missed valence responses. | Confirm in behavioral output; exclude/flag the missing-response component or trials per analysis policy. | [session 8](data_collection_notes/subj6/250118_fMRI_OX_NWU_VS_notes_session8.pdf) |

### Session-level QC

| Session | Run/trial | Note-derived issue | Follow-up | Source note |
|---|---|---|---|---|
| 10 | Run 1 | DAQ trouble, possibly involving the upper USB-C port; note says to make sure run 1's acquisition window matches the others. | Verify trigger timing, duration, and lab-chart alignment before retaining run 1. | [session 10](data_collection_notes/subj6/250118_fMRI_OX_NWU_VS_notes_session10.pdf) |

### Deferred odor-set sanity check (not an automatic exclusion)

Per the project instruction, subject 6 sessions 1 and 2 are recorded separately for a later sanity check. The handwritten notes say the odor list/sets were flipped at collection time and should be fixed before analysis. No trial or run is automatically excluded here on that basis.

- [Session 1 note](data_collection_notes/subj6/250118_fMRI_OX_NWU_VS_notes_session1and2.pdf): odor list flipped; "fix this before data analysis."
- [Session 2 note](data_collection_notes/subj6/250118_fMRI_OX_NWU_VS_notes_session1and2.pdf): odor sets were also flipped at collection time.

No explicit exclusion-worthy event was noted for sessions 7 or 9.

## Recommended verification order

1. Remove or replace explicitly aborted/repeated runs (subject 4 sessions 5 and 11; subject 5 sessions 5, 9, and 10; subject 6 old run 6; incomplete runs noted for subjects 3 and 4).
2. Confirm no-odor, weak-odor, and missed-response trials in behavioral/DAQ/lab-chart records.
3. Resolve motion events to exact timestamps and apply one consistent quantitative motion rule.
4. Resolve session/result/file mappings before concatenating data.
5. Perform the separate subject 6 session 1-2 odor-set sanity check.
