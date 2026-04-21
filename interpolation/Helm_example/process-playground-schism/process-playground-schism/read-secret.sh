#!/bin/bash

function exportSecret() {
    RECIPIENT_VARIABLE=${1}
    SECRET_KEY=${2}
    if [ -z "${!SECRET_KEY}" ]
    then
        local __SECRET
        __SECRET=""
        local __MOI_TOOL
        __MOI_TOOL=$(command -v moi) || true
        if [ "$__MOI_TOOL" != "" ]
        then
            echo "Reading secret $SECRET_KEY with the moi tool..." >&2
            __SECRET="$("$__MOI_TOOL" read-secret --name "$SECRET_KEY")" || true
        fi
        if [ "$__SECRET" = "" ]
        then
            echo "${SECRET_KEY} not found" >&2
        else
            export "${SECRET_KEY}"="$__SECRET"
            echo "${SECRET_KEY} cached" >&2
        fi
    fi
    if [ -z "${!SECRET_KEY}" ]
    then
        echo "${RECIPIENT_VARIABLE} not exported" >&2
    else
        echo "${RECIPIENT_VARIABLE} exported" >&2
        export "${RECIPIENT_VARIABLE}"=${!SECRET_KEY}
    fi
}
